from dataclasses import replace
import numpy as np
import pytest
from src.model import Parameters, flight_model
from src.simulation import simulate
from src.rk4 import rk4_step

Y0 = np.array([30000., 2112.14, np.deg2rad(-2.), 300., 1400.])


@pytest.mark.parametrize('state,expected', [
    ([30000., 2112.14, np.deg2rad(-2.), 300., 1400.],
     [-73.71262296522039, -2.006120178267873, -0.00015681280689557736, 309.09488520177536, 0.]),
    ([20000., 1800., np.deg2rad(-4.), 2200., 1400.],
     [-125.56165273942554, -6.427442492145637, 0.010391643914155193, 234.05040846215658, 0.]),
])
def test_independent_derivative_benchmarks(state, expected):
    # Reference numbers separately evaluated with Python math scalar arithmetic,
    # without calling the project's density, gravity, or heat-flux helpers.
    np.testing.assert_allclose(flight_model(0., np.array(state)), expected,
                               rtol=1e-12, atol=1e-12)


@pytest.mark.parametrize('dt', [0., -.1, np.nan, np.inf])
def test_invalid_time_step(dt):
    with pytest.raises(ValueError):
        simulate(Y0, dt, 1.)


@pytest.mark.parametrize('end', [-1., np.nan, np.inf])
def test_invalid_horizon(end):
    with pytest.raises(ValueError):
        simulate(Y0, .1, end)


@pytest.mark.parametrize('index,value', [(1, 0.), (4, 0.), (3, -1.), (0, np.nan)])
def test_invalid_initial_state(index, value):
    y = Y0.copy(); y[index] = value
    with pytest.raises(ValueError):
        simulate(y, .1, 1.)


def test_initial_event_and_zero_duration():
    for kwargs in [dict(t_final=1., stop_altitude=31000.), dict(t_final=0.)]:
        t, y = simulate(Y0, .1, **kwargs)
        assert np.array_equal(t, [0.])
        np.testing.assert_array_equal(y[0], Y0)


def test_partial_final_step_and_input_unchanged():
    original = Y0.copy()
    t, y = simulate(Y0, .3, 1.)
    assert t[-1] == 1. and np.all(np.diff(t) > 0)
    np.testing.assert_array_equal(Y0, original)
    np.testing.assert_array_equal(y[:, 4], np.full(len(t), Y0[4]))


def test_altitude_event_located_inside_step():
    t, y = simulate(Y0, 1., 10., stop_altitude=29999.)
    assert 0 < t[-1] < 1
    assert abs(y[-1, 0]-29999.) < 1e-7


def test_rk4_nonautonomous_vector_equation():
    got = rk4_step(lambda t, y: np.array([t, 2*t]), 2., np.array([2., 4.]), .3)
    np.testing.assert_allclose(got, [2.3**2/2, 2.3**2], rtol=1e-14)


def test_radiative_cooling_exact_limit_fourth_order():
    p = replace(Parameters(), rho30=0., g0=0.)
    y0 = Y0.copy(); y0[2] = 0.; y0[3] = 2500.
    errors = []
    for dt in [.0625, .03125, .015625]:
        t, y = simulate(y0, dt, 20., params=p)
        exact = (y0[3]**-3+3*p.epsilon*p.sigma*t/p.CA)**(-1/3)
        errors.append(np.max(np.abs(y[:, 3]-exact)))
        np.testing.assert_allclose(y[:, :3], np.tile(y0[:3], (len(t), 1)), atol=1e-12)
    assert all(3.6 < np.log2(a/b) < 4.4 for a, b in zip(errors, errors[1:]))


def test_drag_exact_limit():
    p = replace(Parameters(), g0=0., CL=0., epsilon=0., heating_coefficient=0.)
    y0 = Y0.copy(); y0[2] = 0.
    t, y = simulate(y0, 1., 100., params=p)
    a = .5*p.rho30*p.CD*p.A/y0[4]
    np.testing.assert_allclose(y[:, 1], y0[1]/(1+a*y0[1]*t), atol=1e-7, rtol=0)


def test_drag_limit_fourth_order():
    p = replace(Parameters(), CD=2., g0=0., CL=0., epsilon=0., heating_coefficient=0.)
    y0 = Y0.copy(); y0[2] = 0.
    errors = []
    for dt in [1., .5, .25]:
        t, y = simulate(y0, dt, 20., params=p)
        a = .5*p.rho30*p.CD*p.A/y0[4]
        errors.append(np.max(np.abs(y[:, 1]-y0[1]/(1+a*y0[1]*t))))
    assert all(3.8 < np.log2(a/b) < 4.2 for a, b in zip(errors, errors[1:]))


def test_mechanical_energy_derivative():
    p = Parameters()
    dh, dv, *_ = flight_model(0., Y0, p)
    g = p.g0*(p.RE/(p.RE+Y0[0]))**2
    energy_rate = Y0[1]*dv + g*dh
    drag_power = -.5*p.rho30*p.CD*p.A*Y0[1]**3/Y0[4]
    np.testing.assert_allclose(energy_rate, drag_power, rtol=1e-13)


def test_emissivity_does_not_change_trajectory():
    _, low = simulate(Y0, .1, 20., params=replace(Parameters(), epsilon=.5))
    _, high = simulate(Y0, .1, 20., params=replace(Parameters(), epsilon=1.))
    np.testing.assert_array_equal(low[:, :3], high[:, :3])
    assert high[-1, 3] < low[-1, 3]


@pytest.mark.parametrize('kwargs', [dict(epsilon=1.1), dict(Rn=0.), dict(missing=1.)])
def test_invalid_parameters(kwargs):
    with pytest.raises((ValueError, TypeError)):
        Parameters(**kwargs)
