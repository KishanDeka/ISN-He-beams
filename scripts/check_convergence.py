#!/usr/bin/env python3
"""Check that beam fluxes are converged with respect to the velocity grid.

Compares the default production grid (2 deg, 1 km/s) with a 2x and 4x finer
grid at a few observer longitudes.
"""

import argparse

from _common import add_inflow_args, inflow_from_args
from isnhe import flux_samples


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--obs-lons", type=float, nargs="+",
                    default=[15.0, 45.0, 75.0, 130.0, 300.0])
    ap.add_argument("--beta", type=float, default=1.0e-7)
    add_inflow_args(ap)
    args = ap.parse_args()
    params = inflow_from_args(args)

    grids = [(2.0, 1.0), (1.0, 0.5), (0.5, 0.25)]
    print(f"{'lon':>6} {'grid':>12} {'F_dir':>12} {'F_ind':>12} "
          f"{'dF_dir/F':>10} {'dF_ind/F':>10}")
    for lon in args.obs_lons:
        ref = None
        rows = []
        for dang, dv in grids[::-1]:
            s = flux_samples(lon, params, dang=dang, dv=dv)
            fd, fi = s.total(args.beta, "direct"), s.total(args.beta, "indirect")
            if ref is None:
                ref = (fd, fi)
            rows.append((dang, dv, fd, fi, fd / ref[0] - 1, fi / ref[1] - 1))
        for dang, dv, fd, fi, ed, ei in rows[::-1]:
            print(f"{lon:6.1f} {f'{dang}deg/{dv}':>12} {fd:12.5e} {fi:12.5e} "
                  f"{ed:10.1e} {ei:10.1e}")


if __name__ == "__main__":
    main()
