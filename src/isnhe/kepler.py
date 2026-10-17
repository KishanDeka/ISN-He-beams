"""Back-tracing of hyperbolic Keplerian trajectories to infinity.

ISN He atoms are assumed to move on Keplerian hyperbolae under solar gravity
alone. Solar radiation pressure on He is small, but it is not negligible for
the indirect beam (Kowalska-Leszczynska et al. 2023, ApJ 950, 98); it is
neglected here. For a phase-space point (r, v) we need

1. the velocity the atom had at infinity, ``v_inf`` (to evaluate the
   interstellar distribution function via Liouville's theorem), and
2. the ionization "exposure" accumulated along the way.

Closed-form exposure
--------------------
If the ionization rate falls off as the inverse square of the heliocentric
distance, beta(r) = beta0 (R0/r)^2, and since along a Keplerian orbit the
true anomaly evolves as d(nu)/dt = h / r^2 (h = |r x v|), the exposure is

    integral beta dt = beta0 R0^2 / h * integral d(nu) = beta0 * R0^2 * dnu / h,

where dnu is the angle swept by the atom from infinity to the current point.
The survival probability is therefore exp(-beta0 * X) with X = R0^2 dnu / h,
which depends only on the orbit (not on beta0). This lets us compute the
geometry once and rescale for any ionization rate.

Direct vs indirect
------------------
The incoming asymptote sits at true anomaly -nu_inf, with cos(nu_inf) = -1/e.
The swept angle is dnu = nu + nu_inf. Atoms with dnu < pi reach the point on
the short ("direct") arc; atoms with dnu > pi went around the Sun ("indirect").
"""

import numpy as np

from .constants import GM_SUN, R0_KM


def _cross(a, b):
    return np.cross(a, b, axis=-1)


def _dot(a, b):
    return np.sum(a * b, axis=-1)


def backtrace(r_vec, v_vec, mu=GM_SUN, r0=R0_KM):
    """Map phase-space points to their incoming asymptotic state.

    Parameters
    ----------
    r_vec : array (..., 3)
        Heliocentric positions [km].
    v_vec : array (..., 3)
        Heliocentric velocities [km/s].

    Returns
    -------
    dict with
        ``unbound`` : bool array, True for hyperbolic orbits (only these come
                      from infinity),
        ``v_inf_vec`` : velocity at infinity [km/s] (NaN where bound),
        ``swept`` : swept true-anomaly angle from infinity [rad],
        ``X`` : exposure factor R0^2 * swept / h [s]; survival = exp(-beta0 X),
        ``e`` : eccentricity, ``h`` : specific angular momentum [km^2/s].
    """
    r_vec = np.asarray(r_vec, dtype=float)
    v_vec = np.asarray(v_vec, dtype=float)
    r_vec, v_vec = np.broadcast_arrays(r_vec, v_vec)

    r = np.sqrt(_dot(r_vec, r_vec))
    v2 = _dot(v_vec, v_vec)
    energy = 0.5 * v2 - mu / r

    h_vec = _cross(r_vec, v_vec)
    h = np.sqrt(_dot(h_vec, h_vec))

    # Radial orbits (h -> 0) and bound orbits are excluded; they are a set of
    # measure zero / do not originate at infinity.
    unbound = (energy > 0) & (h > 1e-12 * r * np.sqrt(v2) + 1e-300)

    with np.errstate(invalid="ignore", divide="ignore"):
        e_vec = _cross(v_vec, h_vec) / mu - r_vec / r[..., None]
        e = np.sqrt(_dot(e_vec, e_vec))

        p_hat = e_vec / e[..., None]
        q_hat = _cross(h_vec / h[..., None], p_hat)

        cos_nu = np.clip(_dot(p_hat, r_vec) / r, -1.0, 1.0)
        nu = np.arccos(cos_nu)
        nu = np.where(_dot(r_vec, v_vec) < 0.0, -nu, nu)

        nu_inf = np.arccos(np.clip(-1.0 / e, -1.0, 1.0))
        swept = nu + nu_inf

        # velocity at true anomaly -nu_inf (incoming asymptote)
        coef = (mu / h)[..., None]
        v_inf_vec = coef * (
            np.sin(nu_inf)[..., None] * p_hat + (e - 1.0 / e)[..., None] * q_hat
        )

        X = r0**2 * swept / h

    v_inf_vec = np.where(unbound[..., None], v_inf_vec, np.nan)
    swept = np.where(unbound, swept, np.nan)
    X = np.where(unbound, X, np.nan)

    return {"unbound": unbound, "v_inf_vec": v_inf_vec, "swept": swept, "X": X,
            "e": e, "h": h}


def is_direct(swept):
    """True for direct-beam atoms (swept angle from infinity below pi)."""
    return swept < np.pi
