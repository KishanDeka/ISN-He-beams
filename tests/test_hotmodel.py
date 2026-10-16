"""Sanity checks of the hot-model moments."""

import numpy as np

from isnhe import InflowParams, density_samples, flux_samples
from isnhe.constants import AU_KM


def test_far_field_recovers_lism_density_and_velocity():
    p = InflowParams()
    r = 3000 * AU_KM * p.upwind_unit
    s = density_samples(r, p, dang=2.0, dv=0.5, vmax=50.0)
    n = s.total(0.0)
    assert np.isclose(n, p.density, rtol=3e-3)
    # bulk velocity from the first moment
    dirs = s.grid.unit_vectors()[s.pix]
    u = np.sum(s.w[:, None] * s.speed[:, None] * dirs, axis=0) / s.w.sum()
    assert np.allclose(u, p.bulk_velocity, atol=0.1)


def test_upwind_density_at_1au_matches_literature():
    """Swaczyna et al. (2025) quote ~0.01 cm^-3 upwind at 1 au for solar
    minimum (a ~33% reduction from 0.015 cm^-3). With a total ionization rate
    of ~1e-7 s^-1 the hot model should give a similar value."""
    p = InflowParams()
    s = density_samples(AU_KM * p.upwind_unit, p, dang=2.0, dv=0.5, vmax=90.0)
    n = s.total(1.0e-7)
    assert 0.009 < n < 0.013


def test_density_scales_linearly_and_ionization_reduces():
    p = InflowParams()
    s = flux_samples(130.0, p, dang=2.0, dv=1.0)
    f1 = s.total(1e-7, "direct")
    f2 = s.total(1e-7, "direct", density=2 * p.density)
    assert np.isclose(f2, 2 * f1)
    assert s.total(2e-7, "direct") < f1 < s.total(0.0, "direct")


def test_direct_beam_energy_like_ibex():
    """IBEX-Lo sees the ISN He direct beam at ~130 eV around February
    (observer longitude ~130 deg)."""
    s = flux_samples(130.0, dang=2.0, dv=1.0)
    E = s.mean(s.energy_eV, 1e-7, "direct")
    assert 115 < E < 140


def test_grid_convergence():
    a = flux_samples(45.0, dang=2.0, dv=1.0)
    b = flux_samples(45.0, dang=1.0, dv=0.5)
    for which in ("direct", "indirect"):
        assert np.isclose(a.total(1e-7, which), b.total(1e-7, which), rtol=1e-3)
