from math import isclose

import pytest
import numpy as np
from pandas.util.testing import assert_frame_equal

from ..linear_systems import (TorsionalPendulumSystem, SimplePendulumSystem,
                              MassSpringDamperSystem, BaseExcitationSystem)
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
