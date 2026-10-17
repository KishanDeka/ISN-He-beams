# Interstellar neutral helium beams

A small side project: an analytic hot model of interstellar neutral (ISN) helium at 1 au. It looks at whether the ratio of the indirect to the direct beam (as IMAP-Lo could see it) can tell a higher He ionization rate apart from a higher interstellar He density.

Motivation: New Horizons/SWAP pickup ions suggest the He ionization rate is ~43% higher than models (Antonik et al. 2025). But pickup ions and the direct-beam flux both scale with density × survival, so they can't separate the two on their own. Bzowski et al. (2023) proposed comparing the direct and indirect beams instead, and this repo is my attempt at a quick implementation of that idea.

## Idea

Atoms move on Keplerian hyperbolae (no radiation pressure). For β(r) = β₀ (1 au / r)², the exposure along an orbit is β₀ X with X = (1 au)² Δν / h. Indirect atoms wrap around the Sun, so they pick up more exposure, and

$$\ln\frac{F_\mathrm{ind}}{F_\mathrm{dir}} = \mathrm{const} - \beta_0\,\big(\langle X\rangle_\mathrm{ind}-\langle X\rangle_\mathrm{dir}\big)$$

The density (and the detector efficiency) cancels in the ratio.

![trajectories](figures/fig1_trajectories.png)

## Some results

With β₀ = 1×10⁻⁷ s⁻¹, direct beam in February (Earth at 130°) and indirect beam in late October (Earth at 30°):

- direct beam: ~129 eV, ⟨X⟩ ≈ 0.67×10⁷ s
- indirect beam: ~106 eV, ⟨X⟩ ≈ 1.83×10⁷ s, ~59° from the Sun (needs the IMAP-Lo pivot)
- β₀ × 1.43 → F_ind/F_dir × 0.61
- n × 1.43 → F_ind/F_dir unchanged

![response](figures/fig5_exposure_response.png)
![degeneracy](figures/fig6_degeneracy.png)

Figures (PNG + PDF in `figures/`):

1. trajectories (cold model sketch)
2. density along Earth's orbit
3. sky maps of both beams
4. beam fluxes / energies / angles over the year
5. exposure and response to the two hypotheses
6. density vs ionization degeneracy
7. energy threshold

## Not included (yet)

- radiation pressure
- other near-Sun losses (e.g. elastic scattering)
- secondary He ("Warm Breeze")
- time-dependent / non-r⁻² ionization
- instrument response (collimator, pivot schedule, efficiency vs energy)
- systematics from the inflow parameters and the required counting statistics

## Running it

```bash
pip install -e ".[dev]"
pytest
bash scripts/run_all.sh        # makes data/
jupyter nbconvert --to notebook --execute --inplace notebooks/plot_figures.ipynb
```

Quick example:

```python
from isnhe import InflowParams, flux_samples

p = InflowParams()
s = flux_samples(30.0, p)      # Earth at 30 deg (late October)
for beta0 in (1.0e-7, 1.43e-7):
    Fd, Fi = s.total(beta0, "direct"), s.total(beta0, "indirect")
    print(beta0, Fd, Fi, Fi / Fd)
```

## References

- Antonik, M., et al. 2025, ApJ, 995, 57
- Bzowski, M., Kubiak, M. A., Möbius, E., & Schwadron, N. A. 2023, ApJS, 265, 24
- Kubiak, M. A., Bzowski, M., et al. 2023, arXiv:2309.06931
- McComas, D. J., et al. 2015, ApJS, 220, 22
- Swaczyna, P., et al. 2022, ApJS, 259, 42
- Swaczyna, P., Bzowski, M., & Kubiak, M. A. 2025, ApJS, 277, 17

## License

MIT
