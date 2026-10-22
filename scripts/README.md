# Scripts

These command-line scripts run the `isnhe` model (see [src/README.md](../src/README.md)) and save the results to `data/` as compressed `.npz` files. The notebook [notebooks/plot_figures.ipynb](../notebooks/plot_figures.ipynb) only reads these files and makes plots, so the slow part runs once here.

| Script | Output | Used for |
|---|---|---|
| [`run_all.sh`](run_all.sh) | everything below | regenerating all data in one go |
| [`compute_annual_beams.py`](compute_annual_beams.py) | `data/annual_beams.npz` | Figs 4, 5, 6, 7 and the summary numbers |
| [`compute_skymaps.py`](compute_skymaps.py) | `data/skymaps.npz` | Fig 3 |
| [`compute_density_profile.py`](compute_density_profile.py) | `data/density_profile.npz` | Fig 2 |
| [`check_convergence.py`](check_convergence.py) | `data/convergence.txt` (via `tee`) | checking that the grid is fine enough |
| [`_common.py`](_common.py) | | shared helpers, not run directly |

Fig 1 (trajectories) is computed directly in the notebook from the cold-model formula, so it needs no data file.

Every Python script takes `--help`. They all accept the same inflow-parameter options, so you can rerun everything for different interstellar conditions, for example:

```bash
python scripts/compute_annual_beams.py --temperature 8000 --speed 25.9 --out data/annual_hot.npz
```

---

## `run_all.sh`

```bash
bash scripts/run_all.sh
```

Runs the four scripts below in order with their default settings. It moves to the repo root first, so it works from any directory. `set -euo pipefail` stops it at the first error. It takes about 15 minutes on a laptop, and almost all of that is `compute_annual_beams.py`.

---

## `_common.py`: shared helpers

- **Path setup.** `REPO` is the repo root, and `src/` is added to `sys.path`, so the scripts run even without `pip install -e .`. `DATA` points to `data/`.
- **`add_inflow_args(parser)`** adds `--density`, `--speed`, `--lon`, `--lat` and `--temperature`, with defaults taken from `InflowParams()`: the IBEX consensus values of McComas et al. (2015), with $n = 0.015$ cm⁻³. These five numbers fully describe the interstellar He far from the Sun, a drifting Maxwellian.
- **`add_grid_args(parser, dang, dv, vmax)`** adds the velocity-space grid options:
  - `--dang`: sky pixel size in degrees, i.e. the angular step of velocity directions;
  - `--dv`: speed step in km/s;
  - `--vmax`: largest speed included, in km/s.

  Each script passes its own defaults: coarser for long runs, finer for sky maps.
- **`inflow_from_args` / `grid_from_args`** turn the parsed options into an `InflowParams` object and a dict of grid keywords for `flux_samples` / `density_samples`.
- **`save_npz(path, meta, **arrays)`** creates the folder if needed and saves the arrays with `np.savez_compressed`. It also stores a JSON string `meta` with the inflow parameters, the grid and the units, so every data file records how it was made. The notebook reads `meta` back to rebuild `InflowParams` and the sky grid.

---

## `compute_annual_beams.py`: the beams over a whole year

**Physics.** As Earth goes around the Sun, the direction and energy at which ISN He arrives change a lot, and the indirect beam is only visible in some seasons. The script puts an observer at 1 au at every ecliptic longitude (every 2.5° by default, 144 positions). For each position it computes the sky-integrated flux of the direct and indirect beams as a function of the ionization rate $\beta_0$.

**How.** It calls `isnhe.annual_beam_fluxes`, which for each longitude:

1. builds `flux_samples` once, the expensive step that traces every velocity cell back to infinity;
2. evaluates $F(\beta_0) = \sum_i w_i e^{-\beta_0 X_i}$ for all $\beta_0$ at once, separately for direct ($\Delta\nu < 180°$) and indirect ($\Delta\nu > 180°$) atoms, and for each energy threshold;
3. computes the mean energy and sky centroid of each beam at `--beta-ref`.

Since $\beta_0$ only enters through $e^{-\beta_0X}$, the whole $\beta_0$ table costs almost nothing extra.

**Options.**

| Option | Default | Meaning |
|---|---|---|
| `--lon-step` | 2.5 | observer longitude step (deg) |
| `--beta-max`, `--nbeta` | 3×10⁻⁷, 61 | $\beta_0$ grid from 0 to `beta-max` s⁻¹ (steps of 5×10⁻⁹) |
| `--beta-ref` | 1×10⁻⁷ | $\beta_0$ used for energies and centroids |
| `--emins` | 0 20 40 | energy thresholds in the observer frame (eV), a crude stand-in for the detector's low-energy efficiency cut-off |
| grid | 2°, 1 km/s, 120 km/s | coarser than the sky maps; `check_convergence.py` shows this is enough |

**Output `data/annual_beams.npz`.**

| Key | Shape | Content |
|---|---|---|
| `obs_lons` | (144,) | observer ecliptic longitude (deg) |
| `betas` | (61,) | ionization rate at 1 au (s⁻¹) |
| `emins` | (3,) | energy thresholds (eV) |
| `F_dir`, `F_ind` | (3, 144, 61) | sky-integrated beam flux (cm⁻² s⁻¹) |
| `E_dir`, `E_ind` | (144,) | flux-weighted mean energy in the observer frame (eV) |
| `cen_dir`, `cen_ind` | (144, 3) | beam centroid lon, lat and RMS width (deg) |
| `beta_ref`, `meta` | | reference $\beta_0$, and the run settings as JSON |

The notebook interpolates these tables to any $\beta_0$ in $\ln F$ (`log_interp_beta`) and gets the mean exposure $\langle X\rangle = -d\ln F/d\beta_0$ from them (`effective_exposure`).

---

## `compute_skymaps.py`: where the beams appear on the sky

**Physics.** A detector only sees the directions it points at, so we need the flux per unit solid angle across the sky, not only the total. The script maps the directional flux (cm⁻² s⁻¹ sr⁻¹) over ecliptic longitude and latitude of arrival direction, in the observer's moving frame. Earth's 29.78 km/s orbital motion shifts the beams on the sky (an aberration effect), and it is included.

**How.** For each observer longitude it builds `flux_samples` on a finer grid (1° pixels, 0.5 km/s) and calls `skymap(beta0, "direct" / "indirect")`. That bins the weighted samples by arrival pixel and divides by the pixel solid angle.

**Options.** `--obs-lons` (default 30, 45 and 130°, i.e. late October, early November and early February), `--betas` (default 0.7×10⁻⁷ and 1.0×10⁻⁷ s⁻¹), and `--emin` (energy threshold, default 0).

**Output `data/skymaps.npz`.** `obs_lons`, `betas`, `lon` and `lat` (pixel centres), and `direct` and `indirect` with shape `(n_obs_lon, n_beta, n_lat, n_lon)`, i.e. (3, 2, 180, 360). The notebook draws the total in grey and contours of each beam, then overlays scan circles from `isnhe.scan_circle`: the IBEX-Lo 90° circle, and the IMAP-Lo pivot angle that passes through the indirect beam.

---

## `compute_density_profile.py`: He density along Earth's orbit

**Physics.** This is mostly a check that the model behaves sensibly. Solar gravity bends the incoming atoms' paths so they meet behind the Sun, forming a **focusing cone** downwind (Earth passes through it in early December), where the density is several times higher. Ionization removes atoms, especially those that pass close to the Sun, so a larger $\beta_0$ lowers the density and shrinks the cone. Upwind at 1 au the model should give about 0.01 cm⁻³ for solar minimum (Swaczyna et al. 2025).

**How.** For each point on a 1 au circle in the ecliptic (every 2.5°), it builds `density_samples` in the Sun's rest frame. Density is a frame-independent number, so the observer's motion doesn't matter here. It then evaluates the total, direct and indirect density for each $\beta_0$.

**Options.** `--lon-step` (2.5°), `--betas` (0, 0.7, 1.0 and 1.43 ×10⁻⁷ s⁻¹: no ionization, low, model, and the 43% higher rate from New Horizons/SWAP). The grid defaults to 2° and 1 km/s, with `vmax` of 100 km/s, because heliocentric speeds at 1 au are lower than the speeds relative to a moving observer.

**Output `data/density_profile.npz`.** `lons` (deg), `betas` (s⁻¹), and `n_all`, `n_dir` and `n_ind` with shape `(n_lon, n_beta)` in cm⁻³.

---

## `check_convergence.py`: is the grid fine enough?

**Why.** The integrals are done with the midpoint rule on a grid of directions × speeds. If the grid is too coarse, the narrow beams (only a few degrees wide) are poorly sampled and the fluxes are wrong.

**How.** At a few observer longitudes (15, 45, 75, 130 and 300° by default) it computes $F_\mathrm{dir}$ and $F_\mathrm{ind}$ at the production grid (2°, 1 km/s) and at 2× and 4× finer grids. It prints the fractional differences relative to the finest grid.

**Result** (in `data/convergence.txt`): the production grid differs from the finest by about 5×10⁻⁵ in most seasons and by 2–4×10⁻⁴ at 75°, where the two beams overlap and the flux varies most sharply. Both are far below anything that matters for the beam ratio.

**Options.** `--obs-lons`, `--beta`, and the inflow options.

---

## Output data at a glance

```
data/
├── annual_beams.npz     # F_dir/F_ind(E_min, observer lon, beta0), energies, centroids
├── skymaps.npz          # directional flux maps for a few seasons
├── density_profile.npz  # He density along Earth's orbit
└── convergence.txt      # grid convergence table
```

All `.npz` files can be opened with `np.load(path)`, and `json.loads(str(d["meta"]))` gives the settings they were made with.
