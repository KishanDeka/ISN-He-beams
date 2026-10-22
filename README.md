# Interstellar neutral helium beams

A small side project: an analytic hot model of interstellar neutral (ISN) helium at $1\,\mathrm{au}$. It looks at whether the ratio of the indirect to the direct beam (as IMAP-Lo could see it) can tell a higher He ionization rate apart from a higher interstellar He density.

Motivation: New Horizons/SWAP pickup ions suggest the He ionization rate is $\sim 43\%$ higher than models (Antonik et al. 2025). But pickup ions and the direct-beam flux both scale with density $\times$ survival, so they can't separate the two on their own. Bzowski et al. (2023) proposed comparing the direct and indirect beams instead, and this repo is my attempt at a quick implementation of that idea.

For more details on the physics and the figures, see the [report (PDF)](report/report.pdf).

## Motivation

Atoms move on Keplerian hyperbolae under solar gravity (radiation pressure is ignored for now). Two orbits reach every point near the Sun:

- the **direct** orbit, which sweeps less than $180^{\circ}$ from infinity
- the **indirect** orbit, which wraps around the Sun first ($\Delta\nu > 180^{\circ}$)

If the ionization rate falls off as $\beta(r) = \beta_0 (1\,\mathrm{au}/r)^2$, then since $d\nu/dt = h/r^2$ along the orbit, the exposure has a closed form:

$$\int \beta\, dt = \beta_0 \frac{(1\,\mathrm{au})^2 \Delta\nu}{h} \equiv \beta_0 X$$

where $\Delta \nu$ is the angle swept since infinity and $h$ is the angular momentum. Indirect atoms end up with 2 to 4 times more exposure than direct ones, so

$$\ln\frac{F_\mathrm{ind}}{F_\mathrm{dir}} = \mathrm{const} - \beta_0\,\big(\langle X\rangle_\mathrm{ind}-\langle X\rangle_\mathrm{dir}\big)$$

The density cancels in the ratio, and so does the absolute detector efficiency (as long as it's the same for both beams, which it isn't quite, since the beams arrive at different energies). The closed form is checked against numerical orbit integration in the tests.

![trajectories](figures/fig1_trajectories.png)

## How the model works

- The gas far away is a drifting Maxwellian with the IBEX consensus parameters (McComas et al. 2015): $v = 25.4\,\mathrm{km\,s^{-1}}$, $\lambda = 255.7^{\circ}$, $\beta = 5.1^{\circ}$, $T = 7500\,\mathrm{K}$, $n = 0.015\,\mathrm{cm^{-3}}$.
- At any point near the Sun the distribution function comes from Liouville's theorem, $f(\mathbf{r}, \mathbf{v}) = f_\mathrm{LISM}(\mathbf{v}_\infty)\, e^{-\beta_0 X}$, where $\mathbf{v}_\infty$ and $X$ come from tracing the hyperbola back to infinity analytically (`src/isnhe/kepler.py`).
- Densities and fluxes are integrals over a grid of directions $\times$ speeds. Fluxes are in the frame of an observer on a circular $1\,\mathrm{au}$ orbit (moving at $29.78\,\mathrm{km\,s^{-1}}$).
- Each grid cell stores its $X$, so changing $\beta_0$, the density or an energy threshold doesn't need a recompute.
- The mean exposure of a beam is $\langle X\rangle = -d\ln F/d\beta_0$.

## Some results

With $\beta_0 = 1\times10^{-7}\,\mathrm{s^{-1}}$, direct beam in February (Earth at $130^{\circ}$) and indirect beam in late October (Earth at $30^{\circ}$):

- direct beam: $\sim 129\,\mathrm{eV}$, $\langle X\rangle \approx 0.67\times10^{7}\,\mathrm{s}$
- indirect beam: $\sim 106\,\mathrm{eV}$, $\langle X\rangle \approx 1.83\times10^{7}\,\mathrm{s}$, $\sim 59^{\circ}$ from the Sun (needs the IMAP-Lo pivot)
- H1, $\beta_0 \times 1.43$: $F_\mathrm{dir} \times 0.75$, $F_\mathrm{ind} \times 0.46$, so $F_\mathrm{ind}/F_\mathrm{dir} \times 0.61$
- H2, $n \times 1.43$: both fluxes $\times 1.43$, so $F_\mathrm{ind}/F_\mathrm{dir}$ unchanged

So a $43\%$ higher ionization rate drops the ratio by $\sim 40\%$, while more helium doesn't change it at all. One catch is that the two beams arrive $\sim 20\%$ apart in energy ($129$ vs $106\,\mathrm{eV}$), so an energy-dependent efficiency won't fully cancel. Bzowski et al. (2023) suggest picking direct/indirect pairs at matched energies instead.

![response](figures/fig5_exposure_response.png)
![degeneracy](figures/fig6_degeneracy.png)

### Figures

PNG and PDF versions are in `figures/`. Full descriptions are in the [report](report/report.pdf).

1. **Trajectories:** cold-gas ($T = 0$) sketch of the direct and indirect orbits that hit Earth in October and February.
2. **Density at $1\,\mathrm{au}$:** He density along Earth's orbit for a few $\beta_0$ values, split into direct and indirect atoms. A sanity check: the upwind value is $\sim 0.011\,\mathrm{cm^{-3}}$ (Swaczyna et al. 2025 use $\sim 0.01$), and the downwind peak is the gravitational focusing cone.
3. **Sky maps:** where the two beams show up on the sky in October and February, with the IBEX-Lo $90^{\circ}$ scan circle and the pivot angle needed to hit the indirect beam. In February the indirect beam comes from too close to the Sun to see.
4. **Over the year:** beam fluxes, mean energies and angle from the Sun vs Earth's longitude. Grey bands mark where the two beams overlap.
5. **Exposure and response:** $\langle X\rangle$ for each beam, and how $F_\mathrm{dir}$, $F_\mathrm{ind}$ and their ratio respond to H1 and H2.
6. **Degeneracy:** fit of $(n\,\varepsilon,\ \beta_0)$ to noise-free fake data from H1. The direct beam alone gives a degenerate band; adding the indirect beam closes the contour.

## Not included (yet)

- radiation pressure: changes the indirect-beam flux by up to $\sim 9\%$ (Kowalska-Leszczy&nacute;ska et al. 2023), and it breaks the closed form for $X$
- other near-Sun losses (e.g. elastic scattering by the solar wind): these would lower the ratio just like a higher $\beta_0$ does
- secondary He (the "Warm Breeze"), which has its own direct and indirect beams
- time-dependent and non-$r^{-2}$ ionization (e.g. electron impact). The two beams are observed months apart, so this matters.
- instrument response (collimator, pivot schedule, efficiency vs energy)
- systematics from the inflow parameters and the required counting statistics

## Report

A longer write-up of the physics and the figures is in [report/report.pdf](report/report.pdf) (source: `report/report.tex`, build with `bash report/run.sh`).

## Running it

```bash
pip install -e ".[dev]"
pytest
bash scripts/run_all.sh        # makes data/
jupyter nbconvert --to notebook --execute --inplace notebooks/plot_figures.ipynb
```

The code is explained, with the physics behind each part, in [src/README.md](src/README.md) (the model) and [scripts/README.md](scripts/README.md) (the data-generation scripts).

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
- Bzowski, M., Kubiak, M. A., M&ouml;bius, E., & Schwadron, N. A. 2023, ApJS, 265, 24
- Kowalska-Leszczy&nacute;ska, I., Kubiak, M. A., & Bzowski, M. 2023, ApJ, 950, 98
- Kubiak, M. A., Bzowski, M., et al. 2023, arXiv:2309.06931
- McComas, D. J., et al. 2015, ApJS, 220, 22
- Swaczyna, P., et al. 2022, ApJS, 259, 42
- Swaczyna, P., Bzowski, M., & Kubiak, M. A. 2025, ApJS, 277, 17

## License

MIT
