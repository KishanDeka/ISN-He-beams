"""Classical hot model of interstellar neutral helium in the inner heliosphere.

The distribution function at a point (r, v) inside the heliosphere follows
from Liouville's theorem:

    f(r, v) = f_LISM(v_inf(r, v)) * exp(-beta0 * X(r, v)),

where v_inf and X come from :func:`isnhe.kepler.backtrace`. Moments (density,
flux) are obtained by integrating f on a grid in velocity space.

Because beta0 only enters through exp(-beta0 X), the expensive part (orbit
back-tracing on the grid) is done once; the result is stored as a set of
weighted samples (:class:`VelocitySamples`) that can be re-weighted for any
ionization rate. The ISN density enters linearly, so it can also be rescaled
freely afterwards.
"""

from dataclasses import dataclass

import numpy as np

from .constants import AU_KM, EV, KM_TO_CM, M_HE, V_EARTH
from .kepler import backtrace, is_direct
from .lism import InflowParams, ecliptic_unit_vector, maxwellian


# --------------------------------------------------------------------------
# grids
# --------------------------------------------------------------------------

@dataclass(frozen=True)
class SkyGrid:
    """Regular longitude/latitude grid of directions (ecliptic frame)."""

    dang: float  # pixel size [deg]

    @property
    def lon_edges(self):
        return np.arange(0.0, 360.0 + 0.5 * self.dang, self.dang)

    @property
    def lat_edges(self):
        return np.arange(-90.0, 90.0 + 0.5 * self.dang, self.dang)

    @property
    def lon(self):
        e = self.lon_edges
        return 0.5 * (e[1:] + e[:-1])

    @property
    def lat(self):
        e = self.lat_edges
        return 0.5 * (e[1:] + e[:-1])

    @property
    def shape(self):
        return (self.lat.size, self.lon.size)

    def pixel_solid_angle(self):
        """Solid angle of each pixel [sr], shape (nlat, nlon)."""
        le = np.deg2rad(self.lat_edges)
        band = np.deg2rad(self.dang) * (np.sin(le[1:]) - np.sin(le[:-1]))
        return np.repeat(band[:, None], self.lon.size, axis=1)

    def unit_vectors(self):
        """Pixel-centre unit vectors flattened in (lat, lon) C-order, (N, 3)."""
        lon2, lat2 = np.meshgrid(self.lon, self.lat)
        return ecliptic_unit_vector(lon2.ravel(), lat2.ravel())


def observer_state(obs_lon_deg, r_au=1.0, v_orbit=V_EARTH):
    """Position [km] and velocity [km/s] of an observer on a circular,
    prograde orbit in the ecliptic at ecliptic longitude ``obs_lon_deg``."""
    lam = np.deg2rad(obs_lon_deg)
    r_vec = r_au * AU_KM * np.array([np.cos(lam), np.sin(lam), 0.0])
    v_vec = v_orbit * np.array([-np.sin(lam), np.cos(lam), 0.0])
    return r_vec, v_vec


# --------------------------------------------------------------------------
# weighted velocity-space samples
# --------------------------------------------------------------------------

@dataclass
class VelocitySamples:
    """Non-negligible grid cells of the velocity-space integral.

    Each sample i carries a weight ``w[i]`` (for density: f d^3v in cm^-3; for
    flux: f v^3 dv dOmega in cm^-2 s^-1) evaluated *without* ionization, the
    exposure factor ``X[i]`` [s], the swept angle, the pixel index (flux only)
    and the speed relative to the observer.
    """

    w: np.ndarray
    X: np.ndarray
    swept: np.ndarray
    speed: np.ndarray
    pix: np.ndarray
    grid: SkyGrid
    density: float

    # -- selections ------------------------------------------------------
    def _mask(self, which="all", emin_eV=0.0):
        m = np.ones(self.w.size, dtype=bool)
        if which == "direct":
            m &= is_direct(self.swept)
        elif which == "indirect":
            m &= ~is_direct(self.swept)
        elif which != "all":
            raise ValueError("which must be 'all', 'direct' or 'indirect'")
        if emin_eV > 0:
            m &= self.energy_eV >= emin_eV
        return m

    @property
    def energy_eV(self):
        """Kinetic energy relative to the observer [eV]."""
        return 0.5 * M_HE * (self.speed * 1e3) ** 2 / EV

    # -- integrals -------------------------------------------------------
    def total(self, beta0, which="all", emin_eV=0.0, density=None):
        """Integrated moment for ionization rate ``beta0`` [s^-1].

        ``beta0`` may be a scalar or 1-D array. ``density`` rescales the ISN
        density at infinity (default: the one used to build the samples).
        """
        m = self._mask(which, emin_eV)
        w, X = self.w[m], self.X[m]
        scale = 1.0 if density is None else density / self.density
        betas = np.atleast_1d(np.asarray(beta0, dtype=float))
        out = np.array([np.sum(w * np.exp(-b * X)) for b in betas]) * scale
        return out if np.ndim(beta0) else out[0]

    def mean(self, quantity, beta0, which="all", emin_eV=0.0):
        """Weighted mean of a per-sample quantity (e.g. ``energy_eV``)."""
        m = self._mask(which, emin_eV)
        wt = self.w[m] * np.exp(-beta0 * self.X[m])
        return np.sum(wt * np.asarray(quantity)[m]) / np.sum(wt)

    def centroid(self, beta0, which="all", emin_eV=0.0):
        """Flux-weighted mean arrival direction of a beam.

        Returns (lon_deg, lat_deg, rms_deg) where rms_deg is the RMS angular
        distance of the flux from the centroid (a measure of beam width).
        """
        m = self._mask(which, emin_eV)
        wt = self.w[m] * np.exp(-beta0 * self.X[m])
        if wt.sum() <= 0:
            return np.nan, np.nan, np.nan
        dirs = self.grid.unit_vectors()[self.pix[m]]
        c = np.sum(wt[:, None] * dirs, axis=0)
        c /= np.linalg.norm(c)
        ang = np.arccos(np.clip(dirs @ c, -1.0, 1.0))
        rms = np.rad2deg(np.sqrt(np.sum(wt * ang**2) / np.sum(wt)))
        lon = np.rad2deg(np.arctan2(c[1], c[0])) % 360.0
        lat = np.rad2deg(np.arcsin(c[2]))
        return lon, lat, rms

    def skymap(self, beta0, which="all", emin_eV=0.0):
        """Directional flux map [cm^-2 s^-1 sr^-1], shape grid.shape."""
        m = self._mask(which, emin_eV)
        wt = self.w[m] * np.exp(-beta0 * self.X[m])
        npix = self.grid.lon.size * self.grid.lat.size
        flat = np.bincount(self.pix[m], weights=wt, minlength=npix)
        return flat.reshape(self.grid.shape) / self.grid.pixel_solid_angle()


def _sample_velocity_space(r_vec, v_origin, sign, power, params, grid, speeds,
                           dv, chunk, keep_x2):
    """Evaluate f on a (direction x speed) grid around v_origin.

    Velocities are v = v_origin + sign * s * d for pixel directions d and
    speeds s. The integrand weight is f * s^power * ds * dOmega.
    """
    dirs = grid.unit_vectors()
    domega = grid.pixel_solid_angle().ravel()
    vth = params.thermal_speed
    u = params.bulk_velocity
    norm = params.density / (np.pi ** 1.5 * vth**3)

    out = {k: [] for k in ("w", "X", "swept", "speed", "pix")}
    s = speeds[None, :, None]
    for i0 in range(0, dirs.shape[0], chunk):
        d = dirs[i0:i0 + chunk]
        v = v_origin[None, None, :] + sign * s * d[:, None, :]
        bt = backtrace(r_vec[None, None, :], v)
        dvinf = bt["v_inf_vec"] - u
        x2 = np.sum(dvinf * dvinf, axis=-1) / vth**2
        good = bt["unbound"] & (x2 < keep_x2)
        if not np.any(good):
            continue
        pix_idx, sp_idx = np.nonzero(good)
        f0 = norm * np.exp(-x2[good])
        sp = speeds[sp_idx]
        out["w"].append(f0 * sp**power * dv * domega[i0 + pix_idx])
        out["X"].append(bt["X"][good])
        out["swept"].append(bt["swept"][good])
        out["speed"].append(sp)
        out["pix"].append(i0 + pix_idx)

    cat = {k: (np.concatenate(v) if v else np.array([])) for k, v in out.items()}
    cat["pix"] = cat["pix"].astype(np.int64)
    return cat


def flux_samples(obs_lon_deg, params: InflowParams = InflowParams(), dang=1.0,
                 dv=0.5, vmax=120.0, r_au=1.0, v_orbit=V_EARTH, chunk=800,
                 keep_x2=30.0):
    """Weighted samples of the ISN He flux seen by an observer.

    The observer sits at heliocentric distance ``r_au`` in the ecliptic at
    longitude ``obs_lon_deg`` and moves on a circular prograde orbit. Pixel
    directions are *arrival* (look) directions in the ecliptic frame; atoms
    arriving from direction d have velocity v = v_obs - s d with s the speed in
    the observer frame. The weight of each sample is the number flux
    f s^3 ds dOmega [cm^-2 s^-1] at unit survival.

    Parameters
    ----------
    dang : float
        Sky pixel size [deg].
    dv : float
        Speed step in the observer frame [km/s].
    vmax : float
        Maximum speed in the observer frame [km/s].
    keep_x2 : float
        Samples with |v_inf - u|^2 / vth^2 above this are dropped
        (exp(-30) ~ 1e-13).
    """
    grid = SkyGrid(dang)
    r_vec, v_obs = observer_state(obs_lon_deg, r_au, v_orbit)
    speeds = (np.arange(int(np.ceil(vmax / dv))) + 0.5) * dv
    cat = _sample_velocity_space(r_vec, v_obs, -1.0, 3, params, grid, speeds,
                                 dv, chunk, keep_x2)
    cat["w"] = cat["w"] * KM_TO_CM  # cm^-3 km/s -> cm^-2 s^-1
    return VelocitySamples(grid=grid, density=params.density, **cat)


def density_samples(r_vec, params: InflowParams = InflowParams(), dang=1.0,
                    dv=0.5, vmax=120.0, chunk=800, keep_x2=30.0):
    """Weighted samples of the ISN He density at heliocentric point ``r_vec``
    [km] (solar inertial frame). Weights are f d^3v [cm^-3]."""
    grid = SkyGrid(dang)
    speeds = (np.arange(int(np.ceil(vmax / dv))) + 0.5) * dv
    cat = _sample_velocity_space(np.asarray(r_vec, float), np.zeros(3), 1.0, 2,
                                 params, grid, speeds, dv, chunk, keep_x2)
    return VelocitySamples(grid=grid, density=params.density, **cat)


def density(r_vec, beta0, params: InflowParams = InflowParams(), **kw):
    """ISN He density [cm^-3] at ``r_vec`` [km]: dict with total, direct and
    indirect contributions."""
    s = density_samples(r_vec, params, **kw)
    return {k: s.total(beta0, k) for k in ("all", "direct", "indirect")}
