#!/usr/bin/env python3
"""Sensitivity of the beam fluxes to the assumed ISN He inflow parameters.

The beam fluxes are recomputed (on the same observer-longitude and beta grids)
for one-at-a-time perturbations of the inflow speed, temperature and
direction. The perturbation sizes are deliberately generous test values, not
literature uncertainties; scale the results to the uncertainties you trust.

Note that the ISN density needs no variant: fluxes are exactly proportional to
it, so it cancels from the indirect/direct ratio.

Example
-------
    python scripts/compute_systematics.py --out data/systematics.npz
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
    ap.add_argument("--out", default=str(DATA / "systematics.npz"))
    ap.add_argument("--lon-step", type=float, default=5.0)
    ap.add_argument("--beta-max", type=float, default=3.0e-7)
    ap.add_argument("--nbeta", type=int, default=61)
    ap.add_argument("--dspeed", type=float, default=0.5, help="km/s")
    ap.add_argument("--dtemp", type=float, default=500.0, help="K")
    ap.add_argument("--ddir", type=float, default=0.5, help="deg (lon and lat)")
    add_inflow_args(ap)
    add_grid_args(ap)
    args = ap.parse_args()

    base = inflow_from_args(args)
    grid = grid_from_args(args)
    variants = {
        "baseline": base,
        "speed+": base.with_(speed=base.speed + args.dspeed),
        "speed-": base.with_(speed=base.speed - args.dspeed),
        "temp+": base.with_(temperature=base.temperature + args.dtemp),
        "temp-": base.with_(temperature=base.temperature - args.dtemp),
        "lon+": base.with_(lon=base.lon + args.ddir),
        "lon-": base.with_(lon=base.lon - args.ddir),
        "lat+": base.with_(lat=base.lat + args.ddir),
        "lat-": base.with_(lat=base.lat - args.ddir),
    }
    lons = np.arange(0.0, 360.0, args.lon_step)
    betas = np.linspace(0.0, args.beta_max, args.nbeta)

    F_dir, F_ind = [], []
    for name, p in variants.items():
        t0 = time.time()
        r = annual_beam_fluxes(lons, betas, p, emins=(0.0,), **grid)
        F_dir.append(r["F_dir"][0])
        F_ind.append(r["F_ind"][0])
        print(f"{name:9s} done in {time.time() - t0:.0f} s", flush=True)

    meta = {"variants": {k: v.to_dict() for k, v in variants.items()},
            "grid": grid, "perturbations": {"dspeed": args.dspeed,
                                            "dtemp": args.dtemp,
                                            "ddir": args.ddir}}
    save_npz(args.out, meta, names=np.array(list(variants)), obs_lons=lons,
             betas=betas, F_dir=np.array(F_dir), F_ind=np.array(F_ind))


if __name__ == "__main__":
    main()
