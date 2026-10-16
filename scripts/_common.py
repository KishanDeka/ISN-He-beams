"""Shared helpers for the command-line scripts."""

import json
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
# Allow running the scripts without `pip install -e .`
sys.path.insert(0, str(REPO / "src"))

from isnhe import InflowParams  # noqa: E402

DATA = REPO / "data"


def add_inflow_args(parser):
    d = InflowParams()
    g = parser.add_argument_group("ISN He inflow parameters at infinity")
    g.add_argument("--density", type=float, default=d.density, help="cm^-3")
    g.add_argument("--speed", type=float, default=d.speed, help="km/s")
    g.add_argument("--lon", type=float, default=d.lon, help="upwind ecliptic longitude, deg")
    g.add_argument("--lat", type=float, default=d.lat, help="upwind ecliptic latitude, deg")
    g.add_argument("--temperature", type=float, default=d.temperature, help="K")
    return parser


def add_grid_args(parser, dang=2.0, dv=1.0, vmax=120.0):
    g = parser.add_argument_group("velocity-space grid")
    g.add_argument("--dang", type=float, default=dang, help="sky pixel size, deg")
    g.add_argument("--dv", type=float, default=dv, help="speed step, km/s")
    g.add_argument("--vmax", type=float, default=vmax, help="max speed, km/s")
    return parser


def inflow_from_args(args):
    return InflowParams(density=args.density, speed=args.speed, lon=args.lon,
                        lat=args.lat, temperature=args.temperature)


def grid_from_args(args):
    return {"dang": args.dang, "dv": args.dv, "vmax": args.vmax}


def save_npz(path, meta, **arrays):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    import numpy as np
    np.savez_compressed(path, meta=json.dumps(meta), **arrays)
    print(f"saved {path}")
