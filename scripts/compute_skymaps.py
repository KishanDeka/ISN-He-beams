#!/usr/bin/env python3
"""Sky maps of the direct and indirect ISN He beams seen from 1 au.

Maps are directional fluxes [cm^-2 s^-1 sr^-1] on an ecliptic lon/lat grid of
arrival directions, for a set of observer longitudes and ionization rates.

Example
-------
    python scripts/compute_skymaps.py --obs-lons 45 130 --betas 0.7e-7 1.0e-7
"""

import argparse

import numpy as np

from _common import (DATA, add_grid_args, add_inflow_args, grid_from_args,
                     inflow_from_args, save_npz)
from isnhe import flux_samples


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", default=str(DATA / "skymaps.npz"))
    ap.add_argument("--obs-lons", type=float, nargs="+", default=[30.0, 45.0, 130.0])
    ap.add_argument("--betas", type=float, nargs="+", default=[0.7e-7, 1.0e-7])
    ap.add_argument("--emin", type=float, default=0.0, help="eV")
    add_inflow_args(ap)
    add_grid_args(ap, dang=1.0, dv=0.5)
    args = ap.parse_args()

    params = inflow_from_args(args)
    grid = grid_from_args(args)
    maps_dir, maps_ind = [], []
    for lon in args.obs_lons:
        print(f"observer longitude {lon} deg")
        s = flux_samples(lon, params, **grid)
        maps_dir.append([s.skymap(b, "direct", args.emin) for b in args.betas])
        maps_ind.append([s.skymap(b, "indirect", args.emin) for b in args.betas])

    meta = {"params": params.to_dict(), "grid": grid, "emin": args.emin,
            "units": "cm^-2 s^-1 sr^-1; axes (obs_lon, beta, lat, lon)"}
    save_npz(args.out, meta, obs_lons=np.array(args.obs_lons),
             betas=np.array(args.betas), lon=s.grid.lon, lat=s.grid.lat,
             direct=np.array(maps_dir), indirect=np.array(maps_ind))


if __name__ == "__main__":
    main()
