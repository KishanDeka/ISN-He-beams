"""Observables built from the hot model: beam fluxes over the year, the
indirect-to-direct beam ratio, scan circles and simple counting statistics."""

import numpy as np
from scipy.interpolate import RegularGridInterpolator

from .hotmodel import flux_samples
from .lism import InflowParams, ecliptic_unit_vector


def angular_distance(lon1, lat1, lon2, lat2):
    """Great-circle distance between directions given in degrees [deg]."""
    a = ecliptic_unit_vector(lon1, lat1)
    b = ecliptic_unit_vector(lon2, lat2)
    return np.rad2deg(np.arccos(np.clip(np.sum(a * b, axis=-1), -1.0, 1.0)))


def annual_beam_fluxes(obs_lons, betas, params: InflowParams = InflowParams(),
                       emins=(0.0,), beta_ref=1.0e-7, progress=False,
                       **grid_kw):
    """Direct and indirect beam fluxes (integrated over the sky) along the
    observer orbit, tabulated on a grid of ionization rates.

    Parameters
    ----------
    obs_lons : array
        Observer ecliptic longitudes [deg].
    betas : array
        Ionization rates at 1 au [s^-1].
    emins : sequence of float
        Energy thresholds in the observer frame [eV]; one table per value.
    beta_ref : float
        Ionization rate used for the beam energies and centroids.
    **grid_kw
        Passed to :func:`isnhe.hotmodel.flux_samples` (dang, dv, vmax, ...).

    Returns
    -------
    dict with
        ``F_dir``, ``F_ind`` : (n_emin, n_lon, n_beta) beam fluxes [cm^-2 s^-1],
        ``E_dir``, ``E_ind`` : (n_lon,) flux-weighted mean energies [eV],
        ``cen_dir``, ``cen_ind`` : (n_lon, 3) centroid lon, lat, rms width [deg],
        plus the input grids.
    """
    obs_lons = np.asarray(obs_lons, float)
    betas = np.asarray(betas, float)
    emins = np.atleast_1d(np.asarray(emins, float))
    shape = (emins.size, obs_lons.size, betas.size)
    F_dir, F_ind = np.zeros(shape), np.zeros(shape)
    E_dir, E_ind = np.zeros(obs_lons.size), np.zeros(obs_lons.size)
    cen_dir, cen_ind = np.zeros((obs_lons.size, 3)), np.zeros((obs_lons.size, 3))
    for i, lon in enumerate(obs_lons):
        s = flux_samples(lon, params, **grid_kw)
        for k, em in enumerate(emins):
            F_dir[k, i] = s.total(betas, "direct", em)
            F_ind[k, i] = s.total(betas, "indirect", em)
        E_dir[i] = s.mean(s.energy_eV, beta_ref, "direct")
        E_ind[i] = s.mean(s.energy_eV, beta_ref, "indirect")
        cen_dir[i] = s.centroid(beta_ref, "direct")
        cen_ind[i] = s.centroid(beta_ref, "indirect")
        if progress:
            print(f"  obs lon {lon:6.1f} deg: F_dir={F_dir[0, i, 0]:.3e} "
                  f"F_ind={F_ind[0, i, 0]:.3e} cm^-2 s^-1 (beta=0)", flush=True)
    return {"obs_lons": obs_lons, "betas": betas, "emins": emins,
            "beta_ref": beta_ref, "F_dir": F_dir, "F_ind": F_ind,
            "E_dir": E_dir, "E_ind": E_ind, "cen_dir": cen_dir,
            "cen_ind": cen_ind}


def log_interp_beta(betas, F, beta0):
    """Interpolate a flux table F(..., beta) along its last axis in log space.

    ln F is a smooth, nearly linear function of beta (log of a sum of
    exponentials), so linear interpolation of ln F is very accurate.
    """
    lnF = np.log(np.clip(F, 1e-300, None))
    flat = lnF.reshape(-1, lnF.shape[-1])
    out = np.array([np.interp(beta0, betas, row) for row in flat])
    return np.exp(out.reshape(lnF.shape[:-1] + np.shape(beta0)))


def effective_exposure(betas, F, beta0):
    """Flux-weighted mean exposure <X> = -d ln F / d beta0 [s]."""
    lnF = np.log(np.clip(F, 1e-300, None))
    d = np.gradient(lnF, betas, axis=-1)
    flat = d.reshape(-1, d.shape[-1])
    out = np.array([np.interp(beta0, betas, row) for row in flat])
    return -out.reshape(d.shape[:-1] + np.shape(beta0))


# --------------------------------------------------------------------------
# counting statistics
# --------------------------------------------------------------------------

def counts_needed(delta_lnR, nsigma=3.0, count_ratio=1.0):
    """Indirect-beam counts needed to detect a change ``delta_lnR`` of
    ln(F_ind/F_dir) at ``nsigma`` significance (statistics only).

    With Poisson counts N_ind and N_dir = N_ind / count_ratio, the error of
    ln(N_ind/N_dir) is sqrt(1/N_ind + 1/N_dir) = sqrt((1 + count_ratio)/N_ind).
    """
    return (1.0 + count_ratio) * (nsigma / np.abs(delta_lnR)) ** 2


def poisson_deviance(model, observed):
    """Baker & Cousins (1984) Poisson likelihood-ratio statistic, summed."""
    model = np.asarray(model, float)
    observed = np.asarray(observed, float)
    return np.sum(2.0 * (model - observed + observed * np.log(observed / model)),
                  axis=-1)


# --------------------------------------------------------------------------
# scan circles (IBEX-Lo / IMAP-Lo-like geometry)
# --------------------------------------------------------------------------

def scan_circle(skymap, grid, obs_lon_deg, pivot_deg, n_points=720,
                spin_axis_offset_deg=0.0):
    """Sample a sky map along the circle viewed by a spinning detector.

    The spin axis points to the Sun (optionally rotated by
    ``spin_axis_offset_deg`` along the ecliptic), and the boresight makes the
    angle ``pivot_deg`` with it (90 deg corresponds to IBEX-Lo; IMAP-Lo's
    pivot platform can change this angle).

    Returns (spin_phase_deg, look_lon_deg, look_lat_deg, flux). Spin phase 0
    is the northernmost point of the circle.
    """
    sun_dir_lon = (obs_lon_deg + 180.0 + spin_axis_offset_deg) % 360.0
    a = ecliptic_unit_vector(sun_dir_lon, 0.0)
    north = np.array([0.0, 0.0, 1.0])
    e1 = north - np.dot(north, a) * a
    e1 /= np.linalg.norm(e1)
    e2 = np.cross(a, e1)
    phase = np.linspace(0.0, 360.0, n_points, endpoint=False)
    ph = np.deg2rad(phase)[:, None]
    p = np.deg2rad(pivot_deg)
    look = np.cos(p) * a + np.sin(p) * (np.cos(ph) * e1 + np.sin(ph) * e2)
    lat = np.rad2deg(np.arcsin(np.clip(look[:, 2], -1, 1)))
    lon = np.rad2deg(np.arctan2(look[:, 1], look[:, 0])) % 360.0

    # periodic in longitude: pad one column on each side
    L = np.concatenate([[grid.lon[-1] - 360.0], grid.lon, [grid.lon[0] + 360.0]])
    M = np.concatenate([skymap[:, -1:], skymap, skymap[:, :1]], axis=1)
    interp = RegularGridInterpolator((grid.lat, L), M, bounds_error=False,
                                     fill_value=None)
    flux = interp(np.column_stack([lat, lon]))
    return phase, lon, lat, flux
