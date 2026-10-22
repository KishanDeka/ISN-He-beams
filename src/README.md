# `isnhe`: the model code

This folder holds the Python package `isnhe`, a small "hot model" of interstellar neutral (ISN) helium near the Sun. Given the state of the gas far from the Sun, it computes the He density and the He flux an observer at $1\,\mathrm{au}$ would see, split into the **direct** and **indirect** beams.

The package has five modules, each building on the previous ones:

| Module | What it does |
|---|---|
| [`constants.py`](isnhe/constants.py) | physical constants and units |
| [`lism.py`](isnhe/lism.py) | the gas far from the Sun (a drifting Maxwellian) |
| [`kepler.py`](isnhe/kepler.py) | traces an orbit back to infinity and computes the ionization exposure |
| [`hotmodel.py`](isnhe/hotmodel.py) | integrates over velocity space to get densities, fluxes and sky maps |
| [`observables.py`](isnhe/observables.py) | quantities over the year, interpolation in $\beta_0$, scan circles, statistics |

[`__init__.py`](isnhe/__init__.py) re-exports the commonly used names, so `from isnhe import flux_samples` works.

The big picture is:

1. Far away, the He gas is a Maxwellian moving towards the Sun.
2. Each atom then moves on a Keplerian hyperbola (gravity only) and can be ionized on the way.
3. By Liouville's theorem, the phase-space density at a point near the Sun equals the density at infinity on the same orbit, times the probability of surviving ionization.
4. Integrating that over velocities gives densities and fluxes.

---

## Units and coordinates

All internal quantities use:

- distances in $\mathrm{km}$, velocities in $\mathrm{km\,s^{-1}}$, time in $\mathrm{s}$
- ionization rates in $\mathrm{s^{-1}}$
- phase-space density $f$ in $\mathrm{cm^{-3}\,(km/s)^{-3}}$, so that $\int f\,d^3v$ comes out in $\mathrm{cm^{-3}}$

Directions are in the **heliocentric ecliptic frame**: $x$ points to ecliptic longitude $0^{\circ}$ (the vernal equinox), $z$ to the north ecliptic pole. A direction with longitude $\lambda$ and latitude $b$ is the unit vector

$$\hat{\mathbf{n}} = (\cos b\cos\lambda,\ \cos b\sin\lambda,\ \sin b).$$

Earth (the observer) at ecliptic longitude $\lambda_E$ sees the Sun at $\lambda_E + 180^{\circ}$.

---

## `constants.py`

| Name | Value | Meaning |
|---|---|---|
| `GM_SUN` | $1.327\times10^{11}\,\mathrm{km^3\,s^{-2}}$ | solar gravitational parameter $\mu = GM_\odot$ |
| `AU_KM` | $1.496\times10^{8}\,\mathrm{km}$ | astronomical unit |
| `KB` | $1.381\times10^{-23}\,\mathrm{J\,K^{-1}}$ | Boltzmann constant |
| `AMU`, `M_HE` | | atomic mass unit, He-4 mass ($4.0026\,\mathrm{amu}$) |
| `EV` | $1.602\times10^{-19}\,\mathrm{J}$ | electronvolt, for converting kinetic energy |
| `V_EARTH` | $29.78\,\mathrm{km\,s^{-1}}$ | mean orbital speed of Earth (circular orbit) |
| `KM_TO_CM` | $10^{5}$ | converts the flux integral from $\mathrm{cm^{-3}\,km\,s^{-1}}$ to $\mathrm{cm^{-2}\,s^{-1}}$ |
| `R0_KM` | $1\,\mathrm{au}$ | reference distance for the ionization rate, $\beta(r) = \beta_0 (R_0/r)^2$ |

---

## `lism.py`: the gas far from the Sun

"LISM" is the local interstellar medium. Far from the Sun (formally at infinity) the He atoms have not yet felt solar gravity or ionization, and their velocities follow a **drifting Maxwellian**.

### `ecliptic_unit_vector(lon_deg, lat_deg)`

Converts ecliptic longitude and latitude (degrees) to unit vectors using the formula above. It works on arrays, so whole sky grids are converted at once.

### `InflowParams`

A frozen dataclass with the five numbers describing the gas at infinity. The defaults are the IBEX consensus values of McComas et al. (2015):

| Field | Default | Meaning |
|---|---|---|
| `density` | $0.015\,\mathrm{cm^{-3}}$ | He number density at infinity, $n_\mathrm{ISN}$ |
| `speed` | $25.4\,\mathrm{km\,s^{-1}}$ | bulk inflow speed $u$ |
| `lon`, `lat` | $255.7^{\circ}$, $5.1^{\circ}$ | ecliptic direction the gas comes **from** (upwind) |
| `temperature` | $7500\,\mathrm{K}$ | gas temperature $T$ |

Derived properties:

- `thermal_speed`: the most probable thermal speed $v_\mathrm{th} = \sqrt{2k_BT/m_\mathrm{He}}$, about $5.6\,\mathrm{km\,s^{-1}}$ for $7500\,\mathrm{K}$. This is the width of the Maxwellian.
- `upwind_unit`: unit vector towards the upwind direction.
- `bulk_velocity`: $\mathbf{u} = -u\,\hat{\mathbf{n}}_\mathrm{upwind}$. The minus sign is there because the gas comes *from* upwind, so it moves *towards* downwind.
- `with_(**kw)` returns a copy with some fields changed, and `to_dict()` is used to store the parameters in the data files.

### `maxwellian(v_vec, params, density=None)`

The distribution function at infinity:

$$f_\mathrm{LISM}(\mathbf{v}) = \frac{n_\mathrm{ISN}}{\pi^{3/2}v_\mathrm{th}^3}\exp\!\left(-\frac{|\mathbf{v}-\mathbf{u}|^2}{v_\mathrm{th}^2}\right).$$

The prefactor normalizes it so that $\int f\,d^3v = n_\mathrm{ISN}$. The optional `density` argument overrides $n_\mathrm{ISN}$.

---

## `kepler.py`: orbits and ionization exposure

This is the physics core. For any phase-space point $(\mathbf{r}, \mathbf{v})$ near the Sun, it answers two questions:

1. What velocity $\mathbf{v}_\infty$ did this atom have far from the Sun? This is needed to evaluate $f_\mathrm{LISM}$.
2. How much ionization has it been exposed to on the way in?

Radiation pressure on He is neglected, so the only force is solar gravity and every atom from infinity moves on a **hyperbola** with the Sun at one focus.

### `backtrace(r_vec, v_vec, mu=GM_SUN, r0=R0_KM)`

It is fully vectorized: `r_vec` and `v_vec` can be arrays of shape `(..., 3)`. Step by step:

**1. Energy and angular momentum.**

$$E = \tfrac12 v^2 - \frac{\mu}{r}, \qquad \mathbf{h} = \mathbf{r}\times\mathbf{v}, \quad h = |\mathbf{h}|.$$

Only orbits with $E > 0$ (hyperbolic) come from infinity. Bound orbits ($E<0$) and purely radial ones ($h \to 0$) are flagged `unbound = False`, and all their outputs are set to NaN. Radial orbits are a set of measure zero, so dropping them doesn't change any integral.

**2. Eccentricity vector.**

$$\mathbf{e} = \frac{\mathbf{v}\times\mathbf{h}}{\mu} - \frac{\mathbf{r}}{r}, \qquad e = |\mathbf{e}| > 1 \text{ for a hyperbola}.$$

$\mathbf{e}$ points to perihelion. Together with $\mathbf{h}$ it defines a frame in the orbital plane:

- $\hat{\mathbf{p}} = \mathbf{e}/e$, pointing to perihelion,
- $\hat{\mathbf{q}} = \hat{\mathbf{h}}\times\hat{\mathbf{p}}$, $90^{\circ}$ ahead in the direction of motion.

**3. True anomaly.** The angle from perihelion to the current position is

$$\cos\nu = \frac{\hat{\mathbf{p}}\cdot\mathbf{r}}{r}.$$

`arccos` gives only $0 \le \nu \le \pi$. The sign comes from the radial velocity: if $\mathbf{r}\cdot\mathbf{v} < 0$ the atom is still falling towards perihelion, so $\nu$ is negative.

**4. Asymptote and swept angle.** On a hyperbola $r = p/(1+e\cos\nu)$, so $r\to\infty$ when $\cos\nu = -1/e$. The incoming asymptote is therefore at

$$\nu = -\nu_\infty, \qquad \nu_\infty = \arccos(-1/e).$$

The angle the atom has swept since it was at infinity is

$$\Delta\nu = \nu + \nu_\infty \qquad (\text{stored as } \texttt{swept}).$$

**5. Velocity at infinity.** On a Keplerian orbit the velocity at true anomaly $\nu$ is

$$\mathbf{v}(\nu) = \frac{\mu}{h}\big[-\sin\nu\,\hat{\mathbf{p}} + (e + \cos\nu)\,\hat{\mathbf{q}}\big].$$

Setting $\nu = -\nu_\infty$ (so $\sin\nu = -\sin\nu_\infty$ and $\cos\nu = -1/e$) gives the code's formula:

$$\mathbf{v}_\infty = \frac{\mu}{h}\big[\sin\nu_\infty\,\hat{\mathbf{p}} + (e - 1/e)\,\hat{\mathbf{q}}\big].$$

Its magnitude is $\sqrt{2E}$, as energy conservation requires.

**6. Ionization exposure.** The ionization rate is assumed to fall off like the solar EUV flux, as $\beta(r) = \beta_0(R_0/r)^2$. The fraction of atoms that survive is $\exp(-\int\beta\,dt)$. Kepler's second law gives $d\nu/dt = h/r^2$, so $dt/r^2 = d\nu/h$ and

$$\int_{-\infty}^{t}\beta\,dt' = \beta_0 R_0^2\int\frac{dt'}{r^2} = \beta_0\,\frac{R_0^2\,\Delta\nu}{h} \equiv \beta_0 X.$$

The code stores `X = r0**2 * swept / h`, in seconds. Two things make this useful:

- it is **exact** for an $r^{-2}$ rate on a Keplerian orbit, so no orbit integration is needed;
- $X$ depends only on the orbit, not on $\beta_0$. The geometry is computed once and the survival $e^{-\beta_0 X}$ can be evaluated for any $\beta_0$ later.

The function returns `unbound`, `v_inf_vec`, `swept`, `X`, `e` and `h`.

### `is_direct(swept)`

Returns `swept < pi`. Two hyperbolae with the same speed at infinity pass through every point near the Sun:

- the **direct** one has swept less than $180^{\circ}$ since infinity,
- the **indirect** one has gone around the Sun first, so it has swept more than $180^{\circ}$.

Indirect atoms pass closer to the Sun (small $h$) and have a larger $\Delta\nu$, so their $X$ is several times larger. That difference is the whole basis of the beam-ratio method.

---

## `hotmodel.py`: densities, fluxes and sky maps

By **Liouville's theorem**, the phase-space density is constant along a trajectory, apart from losses. So at any point near the Sun

$$f(\mathbf{r},\mathbf{v}) = f_\mathrm{LISM}\big(\mathbf{v}_\infty(\mathbf{r},\mathbf{v})\big)\,e^{-\beta_0 X(\mathbf{r},\mathbf{v})}.$$

This is the classical "hot model": hot because the gas has a finite temperature, unlike the "cold model" where every atom has the same velocity at infinity. Moments of $f$ (density, flux) are integrals over velocity, done numerically on a grid.

### `SkyGrid(dang)`

A regular longitude $\times$ latitude grid with pixel size `dang` in degrees, covering the whole sky.

- `lon_edges`, `lat_edges`, `lon`, `lat`: pixel edges and centres.
- `pixel_solid_angle()`: the exact solid angle of each pixel, $\Delta\Omega = \Delta\lambda\,(\sin b_2 - \sin b_1)$. Pixels near the poles are smaller, and this accounts for it.
- `unit_vectors()`: pixel-centre directions, flattened in (lat, lon) order.

The same grid is used both for velocity directions in the integral and for sky-map pixels.

### `observer_state(obs_lon_deg, r_au=1.0, v_orbit=V_EARTH)`

Position and velocity of an observer on a circular, prograde orbit in the ecliptic plane:

$$\mathbf{r} = r(\cos\lambda, \sin\lambda, 0), \qquad \mathbf{v}_\mathrm{obs} = v_\mathrm{orb}(-\sin\lambda, \cos\lambda, 0).$$

This stands in for Earth, IBEX or IMAP; IMAP's real position at L1 and the eccentricity of Earth's orbit are ignored. Earth's $29.78\,\mathrm{km\,s^{-1}}$ motion matters a lot: it is comparable to the gas speed, so it shifts where the beams appear on the sky and what energies they arrive with.

### `_sample_velocity_space(...)`, the integration engine

Velocity space is sampled in spherical coordinates around an origin velocity:

$$\mathbf{v} = \mathbf{v}_\mathrm{origin} + \text{sign}\cdot s\,\hat{\mathbf{d}},$$

with directions $\hat{\mathbf{d}}$ from the `SkyGrid` and speeds $s$ at bin midpoints `(i + 0.5) * dv` (midpoint rule). For each $(\hat{\mathbf{d}}, s)$ cell it:

1. calls `backtrace` to get $\mathbf{v}_\infty$, $X$ and $\Delta\nu$;
2. evaluates the Maxwellian factor $\exp(-|\mathbf{v}_\infty-\mathbf{u}|^2/v_\mathrm{th}^2)$;
3. keeps the cell only if the orbit is unbound and $|\mathbf{v}_\infty-\mathbf{u}|^2/v_\mathrm{th}^2 <$ `keep_x2` (default 30, so dropped cells have $f/f_\mathrm{max} < e^{-30} \approx 10^{-13}$). This skips the large empty parts of velocity space;
4. stores the cell's **weight** $w = f_\mathrm{LISM}\,s^{\text{power}}\,\Delta s\,\Delta\Omega$ *without* the ionization factor, plus its $X$, $\Delta\nu$, speed and pixel index.

Directions are processed in chunks (`chunk`) to limit memory use.

### `flux_samples(obs_lon_deg, params, dang, dv, vmax, ...)`

Sets up the integral for the **flux seen by a moving observer**. An atom arriving *from* direction $\hat{\mathbf{d}}$ at speed $s$ relative to the observer has heliocentric velocity

$$\mathbf{v} = \mathbf{v}_\mathrm{obs} - s\,\hat{\mathbf{d}}$$

(sign $-1$ and origin $\mathbf{v}_\mathrm{obs}$). The number of atoms crossing a unit area per unit time is $\int f\,|\mathbf{v}_\mathrm{rel}|\,d^3v_\mathrm{rel}$, and with $d^3v_\mathrm{rel} = s^2\,ds\,d\Omega$ that becomes

$$F = \int\!\!\int f\,s^3\,ds\,d\Omega \qquad (\text{power} = 3).$$

The weights then have units $\mathrm{cm^{-3}\,(km/s)^{-3}} \times \mathrm{(km/s)^{4}} = \mathrm{cm^{-3}\,km\,s^{-1}}$. Multiplying by `KM_TO_CM` converts them to $\mathrm{cm^{-2}\,s^{-1}}$. Here $\hat{\mathbf{d}}$ is the **look direction** (where the detector points), so the samples can be binned straight into sky maps.

Defaults: $1^{\circ}$ pixels, $0.5\,\mathrm{km\,s^{-1}}$ steps, speeds up to $120\,\mathrm{km\,s^{-1}}$. The fastest beam (the February direct beam, about $130\,\mathrm{eV}$) arrives at about $80\,\mathrm{km\,s^{-1}}$ relative to Earth, so $120\,\mathrm{km\,s^{-1}}$ comfortably covers the distribution.

### `density_samples(r_vec, params, ...)` and `density(r_vec, beta0, params)`

These do the same for the **density** at a point, in the Sun's rest frame: origin $\mathbf{0}$, sign $+1$, and $d^3v = s^2\,ds\,d\Omega$ ($\mathrm{power} = 2$). The weights are in $\mathrm{cm^{-3}}$. `density()` is a convenience wrapper that returns the total, direct and indirect densities for one $\beta_0$.

### `VelocitySamples`, the stored result

A dataclass holding, for every kept cell: `w` (weight at unit survival), `X`, `swept`, `speed`, `pix` (pixel index), plus the `grid` and the `density` used. Since the expensive orbit computation is stored, everything below is cheap:

- `_mask(which, emin_eV)` selects `"all"`, `"direct"` or `"indirect"` cells (via `is_direct`), optionally only those above an energy threshold.
- `energy_eV`: kinetic energy in the observer frame, $E = \tfrac12 m_\mathrm{He}s^2$, in eV. This is the energy the detector sees.
- `total(beta0, which, emin_eV, density)`: the integral with ionization switched on,
  $$F(\beta_0) = \sum_i w_i\,e^{-\beta_0 X_i}.$$
  `beta0` can be an array, which is how whole tables over $\beta_0$ are made in one call. Since the density enters $f$ linearly, `density=` simply rescales the result.
- `mean(quantity, beta0, ...)`: a flux-weighted mean, e.g. the mean beam energy $\langle E\rangle = \sum w_i e^{-\beta_0X_i}E_i / \sum w_i e^{-\beta_0X_i}$.
- `centroid(beta0, ...)`: the beam's mean arrival direction (the normalized flux-weighted sum of the unit vectors), and its RMS angular width (the flux-weighted RMS angle from that centroid). The width is used to decide when the two beams overlap on the sky.
- `skymap(beta0, ...)`: bins the weighted samples into sky pixels (`np.bincount`) and divides by each pixel's solid angle, giving the directional flux in $\mathrm{cm^{-2}\,s^{-1}\,sr^{-1}}$.

---

## `observables.py`: things you can compare to data

### `angular_distance(lon1, lat1, lon2, lat2)`

Great-circle angle between two directions, $\arccos(\hat{\mathbf{a}}\cdot\hat{\mathbf{b}})$, in degrees. Used for the angle of a beam from the Sun (which is the boresight angle a Sun-pointing spinner needs) and for the separation between the two beams.

### `annual_beam_fluxes(obs_lons, betas, params, emins, beta_ref, ...)`

Loops the observer around its orbit. For each longitude it builds `flux_samples` once, then fills:

- `F_dir`, `F_ind` with shape `(n_emin, n_lon, n_beta)`: sky-integrated beam fluxes for every energy threshold and every $\beta_0$;
- `E_dir`, `E_ind`: the flux-weighted mean energy of each beam at `beta_ref`;
- `cen_dir`, `cen_ind`: each beam's centroid (lon, lat) and RMS width at `beta_ref`.

This is the main data product, used by `scripts/compute_annual_beams.py`.

### `log_interp_beta(betas, F, beta0)`

Interpolates a flux table to any $\beta_0$ between grid points. It interpolates **$\ln F$** linearly, not $F$. Since $F = \sum w_i e^{-\beta_0X_i}$, $\ln F$ is a smooth, nearly straight function of $\beta_0$ (exactly straight if all cells had the same $X$), so interpolating the log is very accurate.

### `effective_exposure(betas, F, beta0)`

The flux-weighted mean exposure of a beam:

$$\langle X\rangle = -\frac{d\ln F}{d\beta_0} = \frac{\sum w_i X_i e^{-\beta_0X_i}}{\sum w_i e^{-\beta_0X_i}}.$$

It is computed with a finite-difference derivative (`np.gradient`) along the $\beta_0$ table. This quantity drives the method, since

$$\ln\frac{F_\mathrm{ind}}{F_\mathrm{dir}} \approx \text{const} - \beta_0\big(\langle X\rangle_\mathrm{ind} - \langle X\rangle_\mathrm{dir}\big).$$

The density $n_\mathrm{ISN}$ and a common detection efficiency multiply both fluxes, so they drop out of the ratio, while $\beta_0$ does not.

### `counts_needed(delta_lnR, nsigma, count_ratio)`

The number of indirect-beam counts needed to see a change $\Delta\ln R$ in the beam ratio at `nsigma` significance, from Poisson statistics alone. With $N_\mathrm{ind}$ and $N_\mathrm{dir} = N_\mathrm{ind}/\text{count\_ratio}$ counts, the error on $\ln(N_\mathrm{ind}/N_\mathrm{dir})$ is $\sqrt{1/N_\mathrm{ind}+1/N_\mathrm{dir}}$, which gives

$$N_\mathrm{ind} = (1+\text{count\_ratio})\left(\frac{n_\sigma}{|\Delta\ln R|}\right)^2.$$

It isn't used by the current notebook (the required-statistics analysis was removed for now), but it is kept for later.

### `poisson_deviance(model, observed)`

The Baker & Cousins (1984) likelihood-ratio statistic for Poisson counts,

$$\chi^2_\lambda = 2\sum\big[m - n + n\ln(n/m)\big],$$

with $m$ the model counts and $n$ the observed counts. It behaves like $\chi^2$ (so $\Delta\chi^2$ contours have the usual meaning) but stays correct for small counts. The degeneracy figure (Fig. 6) uses it.

### `scan_circle(skymap, grid, obs_lon_deg, pivot_deg, n_points, spin_axis_offset_deg)`

Simulates what a **spinning, Sun-pointing detector** sees:

- the spin axis $\hat{\mathbf{a}}$ points to the Sun (longitude $\lambda_E + 180^{\circ}$, latitude $0$), optionally offset along the ecliptic;
- the boresight makes a fixed angle $p$ (the pivot angle) with the spin axis. It is $90^{\circ}$ for IBEX-Lo, and adjustable for IMAP-Lo, whose detector sits on a pivot platform;
- as the spacecraft spins, the boresight traces a circle of radius $p$ around the Sun:
  $$\hat{\mathbf{l}}(\phi) = \cos p\,\hat{\mathbf{a}} + \sin p\,(\cos\phi\,\hat{\mathbf{e}}_1 + \sin\phi\,\hat{\mathbf{e}}_2).$$

Here $\hat{\mathbf{e}}_1$ is the direction towards ecliptic north projected perpendicular to the spin axis, so spin phase $\phi = 0$ is the northernmost point, and $\hat{\mathbf{e}}_2 = \hat{\mathbf{a}}\times\hat{\mathbf{e}}_1$. The sky map is interpolated along the circle with a `RegularGridInterpolator`; one extra column is copied onto each side so the interpolation wraps around $0^{\circ}/360^{\circ}$. It returns the spin phase, the look longitude and latitude, and the flux along the circle.

The model does not include a collimator response, so this is the flux along the boresight direction, not a count rate.

---

## What the model does not include (yet)

- **Radiation pressure** makes He orbits slightly non-Keplerian and breaks the closed form for $X$.
- **Non-$r^{-2}$ or time-varying ionization**, such as electron impact (which falls off faster than $r^{-2}$) or solar-cycle changes.
- **Secondary He (the "Warm Breeze")**. It can be added as a second `InflowParams` population: build its samples and add the weights.
- **The instrument response**: collimator, energy-dependent efficiency (only a sharp `emin` cut is available) and duty cycle.
- **Real spacecraft geometry**: L1 position and Earth's orbital eccentricity.

See [report/report.pdf](../report/report.pdf) for the full discussion.
