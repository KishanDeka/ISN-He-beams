"""Check the closed-form back-tracing against direct numerical integration."""

import numpy as np
import pytest
from scipy.integrate import solve_ivp

from isnhe.constants import AU_KM, GM_SUN
from isnhe.kepler import backtrace, is_direct


def integrate_back(r0, v0, beta0, r_stop_au=2.0e4):
    """Integrate an orbit backwards in time until r > r_stop_au.

    State: position [au], velocity [km/s], accumulated exposure int beta dt.
    Time is in units of 1e6 s to keep the numbers O(1).
    """
    T = 1.0e6
    mu_au = GM_SUN / AU_KM**3  # au^3 s^-2

    def rhs(t, y):
        r = y[:3]
        rr = np.linalg.norm(r)
        drdt = -y[3:6] / AU_KM * T            # backwards in time
        dvdt = +GM_SUN / (rr * AU_KM) ** 2 * (r / rr) * T  # -(-GM r/r^3)
        dexp = beta0 / rr**2 * T
        return np.concatenate([drdt, dvdt, [dexp]])

    def stop(t, y):
        return np.linalg.norm(y[:3]) - r_stop_au
    stop.terminal = True

    y0 = np.concatenate([np.asarray(r0) / AU_KM, v0, [0.0]])
    sol = solve_ivp(rhs, (0, 1e9), y0, method="DOP853", rtol=1e-11, atol=1e-12,
                    events=stop)
    y = sol.y[:, -1]
    return y[3:6], y[6], y[:3]


CASES = [
    # (position [au], velocity [km/s]) -- direct and indirect-type orbits
    (np.array([1.0, 0.0, 0.0]), np.array([-30.0, 35.0, 3.0])),
    (np.array([0.3, 0.9, 0.05]), np.array([40.0, -20.0, -2.0])),
    (np.array([-0.5, -0.8, 0.0]), np.array([10.0, -45.0, 5.0])),
    (np.array([0.7, -0.7, 0.1]), np.array([-50.0, -5.0, 1.0])),
]


@pytest.mark.parametrize("r_au,v", CASES)
def test_vinf_and_exposure_match_ode(r_au, v):
    beta0 = 1.0e-7
    r = r_au * AU_KM
    bt = backtrace(r, v)
    assert bt["unbound"]

    v_far, exposure, r_far = integrate_back(r, v, beta0)
    # Speed at infinity from energy conservation at the end point
    vinf = bt["v_inf_vec"]
    r_far_km = np.linalg.norm(r_far) * AU_KM
    speed_inf_ode = np.sqrt(np.dot(v_far, v_far) - 2 * GM_SUN / r_far_km)
    assert np.isclose(speed_inf_ode, np.linalg.norm(vinf), rtol=1e-7)
    # Direction converges as O(GM / (r v^2))
    cosang = np.dot(v_far, vinf) / np.linalg.norm(v_far) / np.linalg.norm(vinf)
    assert np.degrees(np.arccos(min(1.0, cosang))) < 0.05

    # Exposure tail beyond r_stop is ~ beta0 R0^2 / (v_inf r_stop): tiny
    assert np.isclose(exposure, beta0 * bt["X"], rtol=1e-3)


def test_energy_and_bound_orbits():
    r = np.array([AU_KM, 0, 0])
    v_escape = np.sqrt(2 * GM_SUN / AU_KM)
    bt = backtrace(r, np.array([0, 0.9 * v_escape, 0]))
    assert not bt["unbound"]
    v = np.array([0, 1.5 * v_escape, 0])
    bt = backtrace(r, v)
    vinf_expected = np.sqrt(np.dot(v, v) - v_escape**2)
    assert np.isclose(np.linalg.norm(bt["v_inf_vec"]), vinf_expected, rtol=1e-12)


def state_from_elements(e, p_km, nu):
    """Position/velocity on a planar orbit (perifocal frame = ecliptic xy)."""
    r = p_km / (1 + e * np.cos(nu))
    pos = r * np.array([np.cos(nu), np.sin(nu), 0.0])
    k = np.sqrt(GM_SUN / p_km)
    vel = k * np.array([-np.sin(nu), e + np.cos(nu), 0.0])
    return pos, vel


@pytest.mark.parametrize("nu_deg", [-120.0, -30.0, 0.0, 30.0, 100.0, 125.0])
def test_swept_angle_and_classification(nu_deg):
    e, p = 1.5, 1.2 * AU_KM
    nu = np.radians(nu_deg)
    pos, vel = state_from_elements(e, p, nu)
    bt = backtrace(pos, vel)
    nu_inf = np.arccos(-1 / e)
    assert np.isclose(bt["swept"], nu + nu_inf, atol=1e-10)
    assert bool(is_direct(bt["swept"])) == (nu + nu_inf < np.pi)
    # h = sqrt(GM p)
    assert np.isclose(bt["h"], np.sqrt(GM_SUN * p), rtol=1e-12)
