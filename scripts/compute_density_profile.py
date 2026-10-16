#!/usr/bin/env python3
"""ISN He density at 1 au around the ecliptic (validation / context).

Shows the gravitational focusing cone downwind of the Sun and how the
ionization rate shapes it, split into direct and indirect contributions.

Example
-------
    python scripts/compute_density_profile.py --out data/density_profile.npz
"""

import argparse

import numpy as np

from _common import (DATA, add_grid_args, add_inflow_args, grid_from_args,
                     inflow_from_args, save_npz)
from isnhe import density_samples
from isnhe.constants import AU_KM
from isnhe.hotmodel import observer_state


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", default=str(DATA / "density_profile.npz"))
    ap.add_argument("--lon-step", type=float, default=2.5)
    ap.add_argument("--betas", type=float, nargs="+",
                    default=[0.0, 0.7e-7, 1.0e-7, 1.43e-7])
    add_inflow_args(ap)
    add_grid_args(ap, vmax=100.0)
    args = ap.parse_args()

    params = inflow_from_args(args)
    grid = grid_from_args(args)
    lons = np.arange(0.0, 360.0, args.lon_step)
    n_all = np.zeros((lons.size, len(args.betas)))
    n_dir = np.zeros_like(n_all)
    n_ind = np.zeros_like(n_all)
    for i, lon in enumerate(lons):
        r_vec, _ = observer_state(lon)
        s = density_samples(r_vec, params, **grid)
        n_all[i] = s.total(args.betas, "all")
        n_dir[i] = s.total(args.betas, "direct")
        n_ind[i] = s.total(args.betas, "indirect")
        print(f"  lon {lon:6.1f}: n = {n_all[i, 0]:.4e} cm^-3 (beta=0)", flush=True)

    meta = {"params": params.to_dict(), "grid": grid, "r_au": 1.0,
            "units": "cm^-3"}
    save_npz(args.out, meta, lons=lons, betas=np.array(args.betas),
             n_all=n_all, n_dir=n_dir, n_ind=n_ind)


if __name__ == "__main__":
    main()
