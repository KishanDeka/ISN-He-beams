"""Interstellar neutral helium in the pristine very local interstellar medium (VLISM).

The gas far from the Sun ("at infinity") is a drifting Maxwellian. Default values
follow the IBEX ISN helium consensus (McComas et al. 2015, ApJS 220, 22) and the
density adopted in Swaczyna et al. (2025, ApJS 277, 17) and Antonik et al.
(2025, ApJ 995, 57).
"""

from dataclasses import dataclass, asdict, replace

import numpy as np

from .constants import KB, M_HE


def ecliptic_unit_vector(lon_deg, lat_deg):
    """Unit vector(s) for ecliptic longitude/latitude in degrees (shape (..., 3))."""
    lon = np.deg2rad(lon_deg)
    lat = np.deg2rad(lat_deg)
    return np.stack(
        [np.cos(lat) * np.cos(lon), np.cos(lat) * np.sin(lon), np.sin(lat)], axis=-1
    )


@dataclass(frozen=True)
class InflowParams:
    """Parameters of ISN helium at infinity.

    Attributes
    ----------
    density : float
        ISN He number density at infinity [cm^-3].
    speed : float
        Inflow speed [km/s].
    lon, lat : float
        Ecliptic longitude/latitude of the *upwind* direction (where the gas
        comes from) [deg].
    temperature : float
        Temperature [K].
    """

    density: float = 0.015
    speed: float = 25.4
    lon: float = 255.7
    lat: float = 5.1
    temperature: float = 7500.0

    @property
    def thermal_speed(self):
        """Most probable speed sqrt(2 kT / m) [km/s]."""
        return np.sqrt(2.0 * KB * self.temperature / M_HE) / 1.0e3

    @property
    def upwind_unit(self):
        return ecliptic_unit_vector(self.lon, self.lat)

    @property
    def bulk_velocity(self):
        """Bulk velocity vector of the gas [km/s]; points *away* from upwind."""
        return -self.speed * self.upwind_unit

    def with_(self, **kwargs):
        """Return a copy with some parameters changed."""
        return replace(self, **kwargs)

    def to_dict(self):
        return asdict(self)


def maxwellian(v_vec, params: InflowParams, density=None):
    """Drifting Maxwellian f(v) at infinity [cm^-3 (km/s)^-3].

    Parameters
    ----------
    v_vec : array (..., 3)
        Velocities [km/s].
    params : InflowParams
    density : float, optional
        Override params.density (useful for unit-density weights).
    """
    n = params.density if density is None else density
    vth = params.thermal_speed
    dv = v_vec - params.bulk_velocity
    x2 = np.sum(dv * dv, axis=-1) / vth**2
    return n / (np.pi ** 1.5 * vth**3) * np.exp(-x2)
