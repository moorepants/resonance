from math import isclose

import pytest
import numpy as np

from ..linear_systems import (_ParametersDict, _MeasurementsDict,
                              _CoordinatesDict, TorsionalPendulumSystem,
                              SimplePendulumSystem)


def test_nonvalid_parameters_key():
    p = _ParametersDict({})
    with pytest.raises(ValueError):
        p['ben stiller'] = 12.0
    with pytest.raises(ValueError):
        p['time'] = 12.0
    p['a'] = 1.0
    if 'a' in p:
        pass
    del p['a']

def test_setting_measurements_item():
    m = _MeasurementsDict({})
    # not allowed to set measurements
    with pytest.raises(ValueError):
        m['a'] = 12.0

def test_setting_coordinates_item():
    m = _CoordinatesDict({})
    with pytest.raises(ValueError):
        m['a '] = 12.0

    with pytest.raises(ValueError):
        m['time'] = 12.0

    m['first_key'] = 12.0
    with pytest.raises(ValueError):
        m['second_key'] = 12.0

    if 'first_key' in m:
        pass
    del m['first_key']

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

    assert isclose(sys.period(), 2.0 * np.pi * np.sqrt(1.0 / 9.81))
