from math import isclose

import pytest
import numpy as np

from ..linear_systems import (TorsionalPendulumSystem, SimplePendulumSystem,
                              MassSpringDamperSystem)
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

    m, c, k = sys._canonical_coefficients()

    assert isclose(m, 2.0)
    assert isclose(c, 1.0)
    assert isclose(k, 8.0)

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

    with pytest.raises(AttributeError):
        sys.plot_configuration()

    with pytest.raises(AttributeError):
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


def test_mass_spring_damper_system():

    sys = MassSpringDamperSystem()

    sys.constants['mass'] = 1.0
    sys.constants['stiffness'] = 100

    # this will have to be baked into SingleDoFLinearSystem, not really
    # editable by the user (for now) because it dependds on all these hidden
    # values: omega, zeta, etc
    def underdamped_particular_solution():
        # does this do the right thing if zeta=0 (also need crit damped and
        # overdamped)
        r = omega / omega_n
        X = (Fo / k) / np.sqrt(2 * zeta * r)**2 + (1 - r**2)
        phi = np.atan2(2 * zeta * r, 1 - r**2)
        return X * np.sin(omega * t - phi)

    # how does the user know where the force/torque is applied to the model?
    # for SDOF is that obvious?
    traj = sys.sinusoidal_forced_response(amplitude, frequency,
                                          final_time, initial_time, sample_rate)

    # should this be analytic output? this is only relevant to underdamped
    ratio_input_output_amp, frequency = \
        sys.frequency_response(
            input_amplitudes,
            input_frequencies,
            ratio_to_nat_freq=False)  # could optionally output r instead of freq

    # can we have the students do a data based frequency response, like sys id
    # or ratio of spectrums? also, what about a simple fft?

    # this would assume solution amplitude * sin(frequency * t + phi)
    traj = sys.sinusoidal_forcing(amplitude, frequency, final_time)

    # a0 / 2 + a1 * cos(w * t) + b2 * sin(w * t)
    traj = sys.sinusoidal_forcing([a1], [b2], frequency, a0=1.0, final_time=5.0)

    # a0 / 2 + a1 * cos(w * t) + a2 * n * cos(w * t) + b1 * sin(w * t) + b2 * n * cos(w * t)
    traj = sys.sinusoidal_forcing([a1, a2], [b1, b2], frequency, a0=1.0, final_time=5.0)
