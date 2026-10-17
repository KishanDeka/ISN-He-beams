#!/usr/bin/env python3
"""Direct and indirect ISN He beam fluxes at 1 au over a full year.

For each observer ecliptic longitude the sky-integrated flux of the direct and
of the indirect beam is tabulated on a grid of ionization rates beta0 (at
1 au), for several energy thresholds in the observer frame. Beam energies and
centroid directions are stored too.

Example
-------
    python scripts/compute_annual_beams.py --lon-step 2.5 \
        --out data/annual_beams.npz
"""

import argparse
import time

import numpy as np

from _common import (DATA, add_grid_args, add_inflow_args, grid_from_args,
                     inflow_from_args, save_npz)
from isnhe import annual_beam_fluxes


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", default=str(DATA / "annual_beams.npz"))
    ap.add_argument("--lon-step", type=float, default=2.5,
                    help="observer longitude step, deg")
    ap.add_argument("--beta-max", type=float, default=3.0e-7, help="s^-1")
    ap.add_argument("--nbeta", type=int, default=61)
    ap.add_argument("--beta-ref", type=float, default=1.0e-7,
                    help="beta0 used for beam energies/centroids, s^-1")
    ap.add_argument("--emins", type=float, nargs="+", default=[0.0, 20.0, 40.0],
                    help="energy thresholds in the observer frame, eV")
    add_inflow_args(ap)
    add_grid_args(ap)
    args = ap.parse_args()

    params = inflow_from_args(args)
    grid = grid_from_args(args)
    lons = np.arange(0.0, 360.0, args.lon_step)
    betas = np.linspace(0.0, args.beta_max, args.nbeta)

    print(f"inflow: {params}")
    print(f"grid: {grid}, {lons.size} observer longitudes, {betas.size} betas")
    t0 = time.time()
    res = annual_beam_fluxes(lons, betas, params, emins=args.emins,
                             beta_ref=args.beta_ref, progress=True, **grid)
    print(f"done in {time.time() - t0:.0f} s")

    meta = {"params": params.to_dict(), "grid": grid, "units": {
        "F": "cm^-2 s^-1 (sky-integrated beam flux)", "E": "eV",
        "betas": "s^-1 at 1 au", "cen": "deg (lon, lat, rms width)"}}
    save_npz(args.out, meta, **res)


if __name__ == "__main__":
    main()
