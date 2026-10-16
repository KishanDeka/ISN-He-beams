# Interstellar neutral helium beams

**Can the ratio of the indirect to the direct beam of interstellar neutral helium, as IMAP-Lo will measure it, distinguish a higher helium ionization rate from a higher interstellar helium density?**

This repository is an independent, lightweight implementation of the direct/indirect-beam method proposed by Bzowski et al. (2023, ApJS 265, 24). It applies the method to one specific question: the enhanced ionization rate of interstellar neutral (ISN) helium inferred from New Horizons/SWAP pickup ions (Antonik et al. 2025). The model is an analytic hot model of ISN He at 1 au with a full test suite. It predicts both beams over the year and quantifies how well their ratio separates the competing explanations.

---

## Motivation

The ISN He ionization rate appears to be higher than solar-proxy models predict:

- Swaczyna et al. (2022, ApJS 259, 42) found that the evolution of ISN He fluxes measured by IBEX-Lo over a solar cycle requires one of three things: a change in detection efficiency, a higher ionization rate than assumed, or an additional signal source.
- Antonik et al. (2025, ApJ 995, 57) determined the ionization rate from New Horizons/SWAP He⁺ pickup ions and found it to be ≈43% higher than the model rates.

Neither observable can separate the ionization rate from the interstellar density on its own:

- **Pickup-ion production** scales as *n*_ISN · β.
- **The direct-beam flux** at IBEX/IMAP-Lo scales as *n*_ISN · ε · exp(−β₀⟨X⟩), where ε is the poorly known absolute detection efficiency.

So "more ionization", "more helium" and "different efficiency" are partly degenerate.

Bzowski et al. (2023) proposed a calibration-free way out: compare fluxes of the direct and indirect beams. Atoms on indirect orbits wrap around the Sun before reaching 1 au, so they accumulate more ionization exposure. Bzowski et al. (2023) also summarize an earlier attempt at this method by Witte et al., which failed because the indirect beam could not be separated from the background. IMAP-Lo's pivot platform makes the indirect beam observable, and Kubiak et al. (2023) designed a pivot-angle schedule that includes calibration-free determination of the He and Ne ionization rates.

This repository asks a narrower question. If the SWAP result is right, how large is the signal in the beam ratio, how many counts are needed, and which systematics compete with it?

## The key relation

In this model ISN He atoms move on Keplerian hyperbolae; radiation pressure is neglected (see [Limitations](#limitations-and-next-steps)). Suppose the ionization rate falls off as β(r) = β₀ (1 au / r)². Along a Keplerian orbit the true anomaly changes as dν/dt = h / r², so the exposure has a closed form:

$$\int \beta\,dt \;=\; \beta_0\,\frac{(1\,\mathrm{au})^2\,\Delta\nu}{h} \;\equiv\; \beta_0 X ,$$

where Δν is the angle swept since infinity and h is the specific angular momentum. Two orbits reach each point near the Sun:

- the **direct** orbit, with Δν < 180°;
- the **indirect** orbit, which wraps around the Sun (Δν > 180°) and accumulates two to four times more exposure.

Therefore

$$\ln\frac{F_\mathrm{ind}}{F_\mathrm{dir}} = \mathrm{const} - \beta_0\,\big(\langle X\rangle_\mathrm{ind}-\langle X\rangle_\mathrm{dir}\big).$$

The density *n*_ISN cancels exactly. The absolute efficiency ε also cancels, as long as the efficiency is the same for both beams. Where the two beams arrive at different energies, an energy-dependent efficiency does not fully cancel; see below.

The closed form is exact for a pure r⁻² rate and Keplerian motion. It is checked against numerical orbit integration in the tests.

![trajectories](figures/fig1_trajectories.png)

## Results

All results use a reference ionization rate β₀ = 1.0×10⁻⁷ s⁻¹ at 1 au. They compare two hypotheses:

- **H1:** the ionization rate is 1.43× higher.
- **H2:** the interstellar density is 1.43× higher.

| Quantity | Value |
|---|---|
| Direct beam, February (Earth at 130°) | Mean energy 129 eV; ⟨X⟩ = 0.67×10⁷ s; seen 90° from the Sun (IBEX-Lo geometry) |
| Indirect beam, late October (Earth at 30°) | Mean energy 106 eV; ⟨X⟩ = 1.83×10⁷ s; seen ≈59° from the Sun (requires the IMAP-Lo pivot) |
| Response to H1 (β₀ × 1.43) | F_dir × 0.75, F_ind × 0.46, **F_ind/F_dir × 0.61** |
| Response to H2 (*n*_ISN × 1.43) | F_dir × 1.43, F_ind × 1.43, **F_ind/F_dir × 1.00** |
| Statistical precision: 10⁴ direct-beam counts plus the matching ≈750 indirect counts | β₀ recovered to about ±2% (1σ, counting statistics only); H2 excluded at Δχ² ≈ 220 |
| Inflow-parameter systematics in ln R (±0.5 km/s, ±500 K, ±0.5° in λ and β; quadrature sum) | ≈0.05, compared with a signal of ≈0.50 |

What these numbers mean:

- **The signal is large.** A 43% higher ionization rate lowers the beam ratio by about 40%, while a higher density leaves it unchanged.
- **Counting statistics are not the bottleneck.** Tens to a few hundred indirect-beam counts suffice, depending on the autumn observing date.
- **The ±2% is a statistical floor, not a forecast.** Bzowski et al. (2023) estimate that the method can determine the absolute ionization rate at 1 au with an accuracy better than the current ≈15%. A realistic error budget is set by the physics and instrument effects listed under Limitations, not by counts.
- **The chosen pair arrives at different energies.** The February direct beam (≈129 eV) and the October indirect beam (≈106 eV) differ by about 20% in energy. An energy-dependent sputtering efficiency therefore cancels only partly. Bzowski et al. (2023) recommend pairs of direct and indirect observations that arrive at the instrument with the same energy, and that is the natural refinement here.
- **The two observations are months apart.** The beams are observed in different seasons. A time-variable ionization rate can bias the result, as Bzowski et al. (2023) also note.

![response](figures/fig5_exposure_response.png)
![degeneracy](figures/fig6_degeneracy.png)
![season](figures/fig7_season_systematics.png)

### All figures

| Figure | Content |
|---|---|
| `fig1_trajectories` | Cold-model sketch of direct and indirect orbits to Earth in October and February |
| `fig2_density_1au` | Validation: ISN He density along Earth's orbit (focusing cone, upwind value) |
| `fig3_skymaps` | Sky maps of both beams with IBEX-Lo and pivoted scan circles |
| `fig4_annual_beams` | Beam fluxes, energies and angle from the Sun over the year |
| `fig5_exposure_response` | ⟨X⟩ of each beam; response of the observables to H1 and H2 |
| `fig6_degeneracy` | (*n*·ε, β₀) constraints: direct beam only vs direct + indirect |
| `fig7_season_systematics` | Signal vs inflow-parameter systematics, and required counts vs observing date |
| `fig8_energy_threshold` | Effect of a detector energy threshold on both beams |

PNG and PDF versions are in `figures/`.

## Method

1. **Interstellar gas at infinity.** The gas is a drifting Maxwellian with the IBEX ISN He consensus parameters (McComas et al. 2015): $v_\infty = 25.4\ \mathrm{km\,s^{-1}}$, $\lambda = 255.7^\circ$, $\beta = 5.1^\circ$, $T = 7500\ \mathrm{K}$ and $n_\mathrm{ISN} = 0.015\ \mathrm{cm^{-3}}$.

2. **Liouville mapping.** At every point the distribution function is

   $$f(\mathbf{r}, \mathbf{v}) = f_\mathrm{LISM}(\mathbf{v}_\infty)\, \exp(-\beta_0 X).$$

   Here $\mathbf{v}_\infty$ and $X$ come from an analytic back-trace of the hyperbola (`isnhe/kepler.py`): the eccentricity vector, the true anomaly, and the velocity on the incoming asymptote at $\nu = -\arccos(-1/e)$.

3. **Moments.** Densities and fluxes are midpoint-rule integrals on a grid of directions $\times$ speeds. The density (in the solar frame) is

   $$n = \int f\, d^3v ,$$

   and the flux (over arrival directions, in the frame of an observer on a circular $1\ \mathrm{au}$ orbit moving at $29.78\ \mathrm{km\,s^{-1}}$) is

   $$F = \int f\, v_\mathrm{rel}^3\, dv_\mathrm{rel}\, d\Omega .$$

   Each grid cell is stored once with its $X$ and swept angle (`VelocitySamples`). Rescaling to any $\beta_0$ (via $\exp(-\beta_0 X)$), any density (linear) or any energy threshold therefore needs no recomputation.

4. **Observables.**
   - Sky-integrated beam fluxes, flux-weighted energies, beam centroids and widths, and scan circles for a Sun-pointing spinner with an arbitrary boresight (pivot) angle.
   - Poisson counting statistics with the Baker & Cousins (1984) likelihood, as used by Antonik et al. (2025).
   - The average exposure of each beam,

     $$\langle X \rangle = -\frac{d \ln F}{d\beta_0} .$$

## Limitations and next steps

The model is scoped as a feasibility study of the beam-ratio signal. A full IMAP-Lo analysis would need the following, roughly in order of importance:

1. **Radiation pressure.** The model assumes Keplerian orbits. Kowalska-Leszczyńska et al. (2023) showed that radiation pressure on He changes IMAP-Lo fluxes by up to ~9% for the indirect beam, likely above the statistical errors. This exceeds the inflow-parameter systematics estimated here (≈5% in ln R), so it must be included before the ratio is interpreted. It breaks the closed form for X; the natural upgrade is numerical orbits with a radiation-pressure factor.
2. **Other near-Sun losses.** The beam ratio measures the *total* extra loss of indirect atoms, not photoionization alone. Any process that removes atoms preferentially near the Sun lowers the ratio in the same way as a higher ionization rate. Elastic scattering of ISN He by solar-wind particles has been proposed as an alternative to a ~40% higher photoionization rate (*IBEX Observations of Elastic Scattering of Interstellar Helium by Solar Wind Particles*, arXiv:2409.11784). Separating it would need its different energy and angular signatures.
3. **Secondary helium (the "Warm Breeze").** Kubiak et al. (2014) identified a slower, warmer secondary He population arriving from a different direction, with its own direct and indirect beams. Its contribution to the autumn indirect beam must be modelled. The code supports this: compute `flux_samples` for a second `InflowParams` and add the weights.
4. **Time-dependent and non-r⁻² ionization.** Indirect atoms spend months near the Sun, and the two beams are observed in different seasons, so solar-cycle variation of the ionization rate matters. Electron-impact ionization falls off faster than r⁻² and is not captured by a single β₀ (Sokół et al. 2019). Both effects break the closed form; a numerical exposure integral along the orbit is the upgrade.
5. **Instrument response and energy matching.** The model uses sky-integrated fluxes. It does not include the IMAP-Lo collimator, the pivot schedule of Kubiak et al. (2023), the spin-axis offset, the energy-dependent sputtering efficiency or the duty cycle. Choosing direct/indirect pairs at matched energies, as recommended by Bzowski et al. (2023), would reduce the efficiency systematic.
6. **Correlated inflow parameters.** Inflow parameters are perturbed one at a time here. In practice, the speed, direction and temperature from IBEX are correlated and are derived from the same direct-beam data. A joint fit of the inflow parameters and β₀ to both beams is the consistent approach.
7. **Geometry.** IMAP's position at L1 (0.99 au) and the eccentricity of Earth's orbit are ignored. Collisional filtration of He in the outer heliosheath is neglected; it is small.

## Installation and reproduction

```bash
git clone <this repo> && cd isn-helium-beams
python -m venv .venv && source .venv/bin/activate
pip install -e ".[dev]"

pytest                                    # ~15 s
bash scripts/run_all.sh                   # regenerates data/ (~20 min)
jupyter nbconvert --to notebook --execute --inplace notebooks/figures.ipynb
```

Every script accepts `--help`. Inflow parameters and grid resolution can be changed from the command line, for example:

```bash
python scripts/compute_annual_beams.py --temperature 8000 --speed 25.9 --out data/annual_hot.npz
```

## Using the package

```python
from isnhe import InflowParams, flux_samples

p = InflowParams()                 # IBEX consensus inflow, n = 0.015 cm^-3
s = flux_samples(30.0, p)          # observer at ecliptic longitude 30 deg (late October)
for beta0 in (1.0e-7, 1.43e-7):
    Fd, Fi = s.total(beta0, "direct"), s.total(beta0, "indirect")
    print(f"beta0={beta0:.2e}  F_dir={Fd:.3e}  F_ind={Fi:.3e}  ratio={Fi/Fd:.4f}")

s.mean(s.energy_eV, 1e-7, "indirect")   # ~106 eV
s.centroid(1e-7, "indirect")            # (lon, lat, rms width) of the beam in deg
s.skymap(1e-7, "indirect")              # cm^-2 s^-1 sr^-1 on a lon/lat grid
```

## References

- Antonik, M., Swaczyna, P., McComas, D. J., Elliott, H. A., & Bzowski, M. 2025, ApJ, 995, 57
- Baker, S., & Cousins, R. D. 1984, NIMPR, 221, 437
- Bzowski, M., Kubiak, M. A., Möbius, E., & Schwadron, N. A. 2023, ApJS, 265, 24
- Kowalska-Leszczyńska, I., Kubiak, M. A., & Bzowski, M. 2023, ApJ, 950, 98
- Kubiak, M. A., Bzowski, M., Sokół, J. M., et al. 2014, ApJS, 213, 29
- Kubiak, M. A., Bzowski, M., et al. 2023, ApJS, *Science Opportunities for IMAP-Lo Observations of Interstellar Neutral Helium, Neon, and Oxygen during a Maximum of Solar Activity* (arXiv:2309.06931)
- McComas, D. J., Bzowski, M., Fuselier, S. A., et al. 2015, ApJS, 220, 22
- Sokół, J. M., Bzowski, M., & Tokumaru, M. 2019, ApJ, 872, 57
- Swaczyna, P., Bzowski, M., & Kubiak, M. A. 2025, ApJS, 277, 17
- Swaczyna, P., Kubiak, M. A., Bzowski, M., et al. 2022, ApJS, 259, 42
- *IBEX Observations of Elastic Scattering of Interstellar Helium by Solar Wind Particles*, arXiv:2409.11784

## License

MIT, see `LICENSE`.
