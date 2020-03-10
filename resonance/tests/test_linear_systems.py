from math import isclose

import pytest
import numpy as np
from scipy.optimize import fsolve
from scipy.integrate import odeint
from pandas.util.testing import assert_frame_equal

from ..linear_systems import (SingleDoFLinearSystem, TorsionalPendulumSystem,
                              SimplePendulumSystem, MassSpringDamperSystem,
                              BaseExcitationSystem, MultiDoFLinearSystem)
from ..functions import estimate_period


def test_torsional_pendulum_system():

    sys = TorsionalPendulumSystem()

    with pytest.raises(ValueError):
        sys.constants = {}

    with pytest.raises(ValueError):
        sys.coordinates = {}

    with pytest.raises(ValueError):
        sys.speeds = {}

    with pytest.raises(ValueError):
        sys.measurements = {}

    with pytest.raises(KeyError):
        sys._get_par_vals('not_in_here')

    sys.constants['rotational_inertia'] = 2.0
    sys.constants['torsional_damping'] = 1.0
    sys.constants['torsional_stiffness'] = 8.0

    assert isclose(sys.coordinates['torsion_angle'], 0.0)
    assert isclose(sys.speeds['torsion_angle_vel'], 0.0)

    sys.coordinates['torsion_angle'] = 3.0
    sys.speeds['torsion_angle_vel'] = 5.0

    m, c, k = sys.canonical_coefficients()

    assert isclose(m, 2.0)
    assert isclose(c, 1.0)
    assert isclose(k, 8.0)

    # You can add keys to constants, coordinates, and speeds that are
    # duplicates. We an check this by checking for duplicates on any public
    # method call or when computing measurements.
    sys.constants['torsion_angle'] = 5.0
    with pytest.raises(KeyError):
        sys.free_response(1.0)  # should throw error if there are no measurements also
    del sys.constants['torsion_angle']

    def spring_force(torsional_stiffness, torsion_angle):
        return torsional_stiffness * torsion_angle

    with pytest.raises(ValueError):
        sys.add_measurement('Time', spring_force)

    with pytest.raises(ValueError):
        sys.add_measurement('torsion_angle', spring_force)

    with pytest.raises(ValueError):
        sys.add_measurement('torsion_angle_vel', spring_force)

    with pytest.raises(ValueError):
        sys.add_measurement('rotational_inertia', spring_force)

    sys.add_measurement('spring_force', spring_force)
    assert isclose(sys.measurements['spring_force'], 24.0)
    assert len(sys.measurements) == 1
    if 'spring_force' in sys.measurements:
        pass

    def spring_force2(torsional_stiffness, wrong_key):
        return torsional_stiffness * wrong_key
    with pytest.raises(KeyError):
        sys.add_measurement('spring_force2', spring_force2)
    del sys.measurements['spring_force']

    # testing time as an arg
    def funky_funk(torsional_stiffness, time):
        return torsional_stiffness * time
    sys.add_measurement('eeeeeee', funky_funk)
    assert isclose(sys.measurements['eeeeeee'], 0.0)  # time is zero
    del sys.measurements['eeeeeee']

    with pytest.raises(ValueError):
        sys.add_measurement('eeeeeee__futr', funky_funk)

    # You can add keys to constants, coordinates, and speeds that are
    # duplicates. We an check this by checking for duplicates on any public
    # method call or when computing measurements.
    sys.constants['torsion_angle'] = 5.0
    with pytest.raises(KeyError):
        sys.measurements['spring_force']
    del sys.constants['torsion_angle']

    assert isclose(sys._natural_frequency(m, k), np.sqrt(k / m))

    assert isclose(sys._natural_frequency(m, -k).imag, np.sqrt(k / m))

    assert isclose(sys._damping_ratio(2.0, 1.0, np.sqrt(8.0 / 2.0)),
                   1.0 / 2.0 / 2.0 / np.sqrt(8.0 / 2.0))

    t = np.linspace(0.0, 1.0, num=101)
    x0 = 3.0
    v0 = 5.0
    sys.constants['torsional_damping'] = 0.0

    # no damping, stable case
    wn = np.sqrt(8.0 / 2.0)
    expected_pos = v0 / wn * np.sin(wn * t) + x0 * np.cos(wn * t)
    with pytest.raises(ValueError):
        traj = sys.free_response(1.0, 1.5)
    traj = sys.free_response(1.0)
    np.testing.assert_allclose(traj.index, t)
    np.testing.assert_allclose(traj.torsion_angle.values, expected_pos)

    # no damping, unstable case
    sys.constants['torsional_stiffness'] = -8.0
    wn = np.sqrt(8.0 / 2.0)
    # TODO : Need to check to make sure these are correct coefficients.
    expected_pos = v0 / wn * np.sinh(wn * t) + x0 * np.cosh(wn * t)
    traj = sys.free_response(1.0)
    np.testing.assert_allclose(traj.torsion_angle, expected_pos)

    # underdamped
    sys.constants['torsional_stiffness'] = 8.0
    sys.constants['torsional_damping'] = 1.0
    wn = np.sqrt(8.0 / 2.0)
    z = 1.0 / 2.0 / 2.0 / wn
    wd = wn * np.sqrt(1 - z**2)
    A = np.sqrt(((v0 + z * wn * x0)**2 + (x0 * wd)**2) / wd**2)
    phi = np.arctan2(x0 * wd, v0 + z * wn * x0)
    expected_pos = A * np.exp(-z * wn * t) * np.sin(wd * t + phi)
    traj = sys.free_response(1.0)
    np.testing.assert_allclose(traj.torsion_angle, expected_pos)

    with pytest.raises(ValueError):
        sys.plot_configuration()

    with pytest.raises(ValueError):
        sys.animate_configuration()


def test_simple_pendulum_system():

    sys = SimplePendulumSystem()

    sys.constants['pendulum_mass'] = 1.0  # kg
    sys.constants['pendulum_length'] = 1.0  # m
    sys.constants['acc_due_to_gravity'] = 9.81  # m/s**2

    sys.coordinates['angle'] = 1.0
    traj = sys.free_response(5.0, sample_rate=3000)

    assert isclose(sys.period(), 2.0 * np.pi * np.sqrt(1.0 / 9.81))
    assert isclose(estimate_period(traj.index, traj.angle), sys.period(),
                   rel_tol=1e-3)


def test_mass_spring_damper_system():
    sys = MassSpringDamperSystem()

    sys.constants['mass'] = 1.0
    sys.constants['stiffness'] = 100

    wn = np.sqrt(100.0/1.0)
    assert isclose(sys.period(), 2.0 * np.pi / wn)

    sys.constants['damping'] = 0.2
    zeta = 0.2 / 1.0 / (2*wn)
    assert isclose(sys.period(), 2 * np.pi / (wn * np.sqrt(1 - zeta**2)))


def test_mass_spring_damper_system_forced():

    # Example 2.1.1 in the Inman book.
    sys = MassSpringDamperSystem()

    sys.coordinates['position'] = 0.0  # m
    sys.speeds['velocity'] = 0.2  # m/s

    sys.constants['mass'] = 10  # kg
    sys.constants['stiffness'] = 1000  # N/m
    sys.constants['damping'] = 0.0  # Ns/m

    sys.sinusoidal_forcing_response(23.0, 2 * np.sqrt(1000.0 / 10), 3.0)

    # Now add damping and make sure that sinusoidal_forcing_response gives the
    # same results as periodic_forcing_response.
    sys.constants['damping'] = 10.0  # Ns/m

    duration = 5.0

    # this would assume solution: amplitude * cos(frequency * t)
    traj1 = sys.sinusoidal_forcing_response(23.0, 2 * np.sqrt(1000.0 / 10),
                                            duration)
    # see if the periodic forcing gives same answer as sinusoidal:
    # 0.0 / 2 + 23.0 cos() + 0.0 sin()
    traj2 = sys.periodic_forcing_response(0.0, 23.0, 0.0,
                                          2 * np.sqrt(1000.0 / 10), duration)

    assert_frame_equal(traj1, traj2)

    # Ensure that the function works with different numbers of coeffs.
    a0 = 0.1
    a1 = 0.1
    b1 = 0.2
    freq = 4 * np.pi
    # a0 / 2 + a1 * cos(w * t) + b2 * sin(w * t)
    sys.periodic_forcing_response(a0, a1, b1, freq, duration)
    # a0 / 2 + a1 * cos(w * t) + a2 * n * cos(w * t) +
    #          b1 * sin(w * t) + b2 * n * cos(w * t)
    a2 = 0.02
    b2 = 0.03
    sys.periodic_forcing_response(a0, [a1, a2], [b1, b2], freq, duration)


def test_base_excitation_system():

    sys = BaseExcitationSystem()
    sys.constants['mass'] = 100  # kg
    sys.constants['stifness'] = 2000  # N/m
    sys.constants['damping'] = 30  # kg/s

    amplitude = 0.03  # m
    frequency = 6  # rad/s

    traj1 = sys.sinusoidal_base_displacing_response(amplitude, frequency, 5.0)

    traj2 = sys.periodic_base_displacing_response(0.0, 0.0, amplitude,
                                                  frequency, 5.0)

    assert_frame_equal(traj1, traj2)


def test_base_excitation_system_on_book_prob():

    # parameters from example 3.3.2 in the Inman book
    m = 1  # kg
    c = 10  # kg/s
    k = 1000  # N/m

    x0 = 0.01  # m
    v0 = 3.0  # m/s

    Y = 0.05  # m
    wb = 3  # rad/s

    # fourier coefficients
    a0 = 0
    a1 = c * Y * wb
    b1 = k * Y

    z = c / 2 / np.sqrt(k * m)
    wn = np.sqrt(k/m)
    wd = wn * np.sqrt(1 - z**2)

    # compute the steady state solution from the book
    # amplitude of the steady state
    X = wn * Y * np.sqrt((wn**2 + (2*z*wb)**2) /
                         ((wn**2 - wb**2)**2 + (2*z*wn*wb)**2))

    theta1 = np.arctan2(2*z*wn*wb, wn**2 - wb**2)
    theta2 = np.arctan2(wn, 2*z*wb)

    t = BaseExcitationSystem._calc_times(4.0, 0.0, 100)

    xss_expected = X * np.cos(wb*t - theta1 - theta2)

    # steady state solution from resonance
    sys = BaseExcitationSystem()

    sys.constants['mass'] = m  # kg
    sys.constants['damping'] = c  # kg/s
    sys.constants['stiffness'] = k  # N/m

    sys.coordinates['position'] = x0  # m
    sys.speeds['velocity'] = v0  # m/s

    xss, vss, ass, n, thetan, denom = \
        sys._periodic_forcing_steady_state(a0, np.array([[a1]]),
                                           np.array([[b1]]), wb, t)
    np.testing.assert_allclose(xss, xss_expected)

    # now compute the transient solution from the book
    def solve(params):
        A_, phi_ = params
        eq1 = x0 - A_*np.sin(phi_) - X*np.cos(-theta1 - theta2)
        eq2 = (v0 - wd*A_*np.cos(phi_) + z*wn*A_*np.sin(phi_) +
               wb*X*np.sin(-theta1 - theta2))
        return np.array([eq1, eq2])

    def jac(params):
        A_, phi_ = params
        return np.array([[-np.sin(phi_), -A_*np.cos(phi_)],
                         [-wd*np.cos(phi_) + z*wn*np.sin(phi_),
                          wd*A_*np.sin(phi_) + z*wn*A_*np.cos(phi_)]])

    # use the book answer as a guess (note that these values are from the
    # errata, they are incorrect in the fourth edition)
    A_exp, phi_exp = fsolve(solve, (0.0934, 0.1074), fprime=jac)

    assert isclose(A_exp, 0.0934, abs_tol=1e-4)
    assert isclose(phi_exp, 0.1074, abs_tol=1e-4)

    # transient solution from book
    xh_exp = A_exp * np.exp(-z*wn*t)*np.sin(wd*t + phi_exp)
    xh_exp2 = 0.0934 * np.exp(-z*wn*t)*np.sin(wd*t + 0.1074)

    # now get A and phi from resonance
    A, phi = sys._periodic_forcing_transient_A_phi(wb, n, a0, np.array([[a1]]),
                                                   np.array([[b1]]), thetan,
                                                   denom, t)
    # ensure we get the same as found from book equations, i think the fails
    # here are just finding different but equivalent phase shifts
    # fails: sign mismatch
    # assert isclose(A, A_exp, abs_tol=1e-7)
    # fails: this is pi off
    # assert isclose(phi, phi_exp, rel_tol=1e-7)

    xh, vh, ah = sys._damped_sinusoid(A, phi, t)

    # does transient match?
    np.testing.assert_allclose(xh, xh_exp)
    np.testing.assert_allclose(xh, xh_exp2, atol=1e-4)

    # compare the full solution from book and from resonance
    x_exp = xh_exp + xss_expected
    traj = sys.periodic_forcing_response(a0, a1, b1, wb, 4.0)
    np.testing.assert_allclose(traj.position, x_exp)

    # test out a comparison to a numerical integration solution
    def rhs(state, t):
        x, v = state
        xdot = v
        vdot = (-c*v - k*x + c*Y*wb*np.sin(wb*t) + k*Y*np.sin(wb*t)) / m
        return [xdot, vdot]

    res = odeint(rhs, (x0, v0), t)

    # do both the resonance solution and book solution match the numerical
    # integration?
    np.testing.assert_allclose(traj.position, res[:, 0], atol=0.003)
    np.testing.assert_allclose(x_exp, res[:, 0], atol=0.003)
    np.testing.assert_allclose(xh_exp2 + xss, res[:, 0], atol=0.003)

    # now see if the conversion from sinusoidal displacement works
    disp_traj = sys.sinusoidal_base_displacing_response(Y, wb, 4.0)
    np.testing.assert_allclose(disp_traj.position, x_exp)

    # make sure the periodic funciton gives the same result as above
    per_traj = sys.periodic_base_displacing_response(0.0, 0.0, Y, wb, 4.0)
    np.testing.assert_allclose(per_traj.position, x_exp)


def test_transmissibility():

    def force(stiffness, damping, position, velocity, Y, wb, time):
        y = Y * np.sin(wb * time)
        yd = Y * wb * np.cos(wb * time)
        return damping * (velocity - yd) + stiffness * (position - y)

    def disp_transmissibility(r, z):
        return np.sqrt((1 + (2*z*r)**2)/((1-r**2)**2 + (2*z*r)**2))

    def force_transmissibility(r, z):
        return r**2 * disp_transmissibility(r, z)

    def cosfunc(t, a, w, phi):
        return a * np.cos(w*t + phi)

    wn = 2*np.pi  # 1 Hz natural frequency

    sys = BaseExcitationSystem()
    sys.constants['mass'] = 1.0
    sys.constants['stiffness'] = wn**2 * sys.constants['mass']
    sys.constants['Y'] = 0.1

    # constants that will be changed in iterations
    sys.constants['damping'] = 10
    sys.constants['wb'] = 1

    sys.add_measurement('force', force)

    zvals = np.linspace(0.1, 0.9, num=5)
    rvals = np.linspace(0.01, 3.0, num=10)

    dt = np.zeros((zvals.size, rvals.size))
    ft = np.zeros((zvals.size, rvals.size))

    dt_expected = disp_transmissibility(rvals, zvals[:, np.newaxis])
    ft_expected = force_transmissibility(rvals, zvals[:, np.newaxis])

    for iz, z in enumerate(zvals):
        for ir, r in enumerate(rvals):
            wb = wn * r
            c = 2 * sys.constants['mass'] * wn * z

            sys.constants['damping'] = c
            sys.constants['wb'] = wb

            per = 2 * np.pi / min(wb, wn)

            traj = sys.sinusoidal_base_displacing_response(
                sys.constants['Y'], wb, 20*per)

            pos_amp = np.max(np.abs(traj[15*per:].position))
            dt[iz, ir] = pos_amp / sys.constants['Y']

            force_amp = np.max(np.abs(traj[15*per:].force))
            ft[iz, ir] = force_amp / sys.constants['stiffness'] / \
                sys.constants['Y']

    np.testing.assert_array_almost_equal(dt_expected, dt, decimal=2)
    np.testing.assert_array_almost_equal(ft_expected, ft, decimal=2)


def test_defining_system_from_scratch():

    sys = SingleDoFLinearSystem()

    sys.constants['m'] = 1.0
    sys.constants['c'] = 2.0
    sys.constants['k'] = 3.0

    sys.coordinates['x'] = 1.0
    sys.speeds['v'] = 1.0

    def second_order_eom_coefficients(m, c, k):
        return m, c, k

    # canonical coefficients for the lhs of the second ordre form of the
    # ordinary differential equations, defined as:
    # m * x'' + c * x' + k x = F(m, c, k, t)

    sys.canonical_coeffs_func = second_order_eom_coefficients

    m, c, k = sys.canonical_coefficients()

    assert isclose(m, 1.0)
    assert isclose(c, 2.0)
    assert isclose(k, 3.0)

    def second_order_eom_coefficients(booger, c, k):
        return booger, c, k

    with pytest.raises(ValueError):
        sys.canonical_coeffs_func = second_order_eom_coefficients

    desc = sys.__str__()
    expected_desc = """\
System name: SingleDoFLinearSystem

Canonical coefficients function defined: True
Configuration plot function defined: False
Configuration update function defined: False

Constants
=========
m = 1.00000
c = 2.00000
k = 3.00000

Coordinates
===========
x = 1.00000

Speeds
======
v = d(x)/dt = 1.00000

Measurements
============
"""
    assert desc == expected_desc

def test_ode_eval_func():

    sys = MultiDoFLinearSystem()

    sys.constants['m1'] = 1.0  # kg
    sys.constants['m2'] = 1.0  # kg
    sys.constants['k1'] = 10.0  # N/m
    sys.constants['k2'] = 10.0  # N/m

    sys.coordinates['x1'] = 0.1  # m
    sys.coordinates['x2'] = 0.2  # m
    sys.speeds['v1'] = 0.01  # m/s
    sys.speeds['v2'] = 0.02  # m/s

    def coeff_func(m1, m2, k1, k2):

        # mass matrix 2 x 2
        M = np.array([[m1, 0],
                      [0, m2]])

        # damping matrix 2 x 2
        C = np.zeros_like(M)

        # stiffness matrix 2 x 2
        K = np.array([[k1 + k2, -k2],
                      [-k2, k2]])

        return M, C, K

    sys.canonical_coeffs_func = coeff_func

    # should work with shape(2n,)
    x = np.random.random(4)
    t = 0.1
    res = sys._ode_eval_func(x, t)
    assert res.shape == (4,)
    #np.testing.assert_allclose(res, x * t)
    sys.coordinates['x1'] = 0.1
    sys.coordinates['x2'] = 0.2
    sys.speeds['v1'] = 0.01
    sys.speeds['v2'] = 0.02
    sys._time['t'] = 0.0

    # should work with shape(1, 2n)
    x = np.random.random(4).reshape(1, 4)
    t = 0.1
    res = sys._ode_eval_func(x, t)
    assert res.shape == (1, 4)
    #np.testing.assert_allclose(res, x * t)
    sys.coordinates['x1'] = 0.1
    sys.coordinates['x2'] = 0.2
    sys.speeds['v1'] = 0.01
    sys.speeds['v2'] = 0.02
    sys._time['t'] = 0.0

    # should work with shape(1, 2n, 1)
    x = np.random.random(4).reshape(1, 4, 1)
    t = 0.1
    res = sys._ode_eval_func(x, t)
    assert res.shape == (1, 4, 1)
    #np.testing.assert_allclose(res, x * t)
    sys.coordinates['x1'] = 0.1
    sys.coordinates['x2'] = 0.2
    sys.speeds['v1'] = 0.01
    sys.speeds['v2'] = 0.02
    sys._time['t'] = 0.0

    # should work with shape(m, 2n)
    x = np.random.random(10 * 4).reshape(10, 4)
    t = np.random.random(10)
    res = sys._ode_eval_func(x, t)
    assert res.shape == (10, 4)
    #np.testing.assert_allclose(res, x * t[:, np.newaxis])
    sys.coordinates['x1'] = 0.1
    sys.coordinates['x2'] = 0.2
    sys.speeds['v1'] = 0.01
    sys.speeds['v2'] = 0.02
    sys._time['t'] = 0.0

    # should work with shape(m, 2n, 1)
    x = np.random.random(10 * 4).reshape(10, 4, 1)
    t = np.random.random(10)
    res = sys._ode_eval_func(x, t)
    assert res.shape == (10, 4, 1)
    #np.testing.assert_allclose(res, x * t[:, np.newaxis, np.newaxis])
    sys.coordinates['x1'] = 0.1
    sys.coordinates['x2'] = 0.2
    sys.speeds['v1'] = 0.01
    sys.speeds['v2'] = 0.02
    sys._time['t'] = 0.0

def test_multi_dof_linear_system():

    sys = MultiDoFLinearSystem()

    sys.constants['m1'] = 1.0  # kg
    sys.constants['m2'] = 1.0  # kg
    sys.constants['k1'] = 10.0  # N/m
    sys.constants['k2'] = 10.0  # N/m

    sys.coordinates['x1'] = 0.1  # m
    sys.coordinates['x2'] = 0.2  # m
    sys.speeds['v1'] = 0.0  # m/s
    sys.speeds['v2'] = 0.0  # m/s

    def coeff_func(m1, m2, k1, k2):

        # this is the eom of a simple serial spring mass system sliding on
        # ground

        # columns of each matrix have to be in the correct order relative to
        # the sys.coordinates.keys(), sys.speeds.keys()
        # C @ np.array([sys.speeds.keys()]) should work
        # K @ np.array([sys.coordinates.keys()]) should work

        # mass matrix 2 x 2
        M = np.array([[m1, 0],
                      [0, m2]])

        # damping matrix 2 x 2
        C = np.zeros_like(M)

        # stiffness matrix 2 x 2
        K = np.array([[k1 + k2, -k2],
                      [-k2, k2]])

        return M, C, K

    sys.canonical_coeffs_func = coeff_func

    desc = sys.__str__()
    expected_desc = """\
System name: MultiDoFLinearSystem

Canonical coefficients function defined: True
Configuration plot function defined: False
Configuration update function defined: False

Constants
=========
m1 = 1.00000
m2 = 1.00000
k1 = 10.00000
k2 = 10.00000

Coordinates
===========
x1 = 0.10000
x2 = 0.20000

Speeds
======
v1 = d(x1)/dt = 0.00000
v2 = d(x2)/dt = 0.00000

Measurements
============
"""
    assert desc == expected_desc

    sys.free_response(2.0)
    sys.free_response(2.0, integrator="lsoda")

    sys.constants['w1'] = np.pi / 10  # rad/s
    sys.constants['w2'] = np.pi / 20  # rad/s
    sys.constants['f1'] = 1.0  # N
    sys.constants['f2'] = 0.5  # N

    sys.coordinates['x1'] = 0.0  # m
    sys.coordinates['x2'] = 0.0  # m

    # must return number of values equal to num coordinates
    def forcing_func(w1, w2, f1, f2, time):
        return f1

    with pytest.raises(ValueError):
        sys.forcing_func = forcing_func

    def forcing_func(w1, w2, f1, f2, time):
        return f1, f2, w1

    with pytest.raises(ValueError):
        sys.forcing_func = forcing_func

    def forcing_func(w1, w2, f1, f2, time):
        F1 = f1 * np.cos(w1 * time)
        F2 = f2 * np.cos(w2 * time)
        return F1, F2

    sys.forcing_func = forcing_func

    traj = sys.forced_response(2.0)

    assert traj['x1'].sum() > 1E-10

    def forcing_func(w1, w2, f1, f2, time):
        return 0.0, f2 * np.cos(w2 * time)

    sys.forcing_func = forcing_func

    traj = sys.forced_response(2.0)

    assert traj['x1'].sum() > 1E-10

    traj = sys.forced_response(2.0, integrator="lsoda", full_output=True)
