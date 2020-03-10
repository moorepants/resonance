import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
import pytest
from pandas.util.testing import assert_frame_equal

from ..nonlinear_systems import (SingleDoFNonLinearSystem,
                                 MultiDoFNonLinearSystem)


def test_single_dof_nonlinear_system():

    class SimplePendulum(SingleDoFNonLinearSystem):
        pass

    sys = SimplePendulum()

    sys.constants['mass'] = 1.0  # kg
    sys.constants['length'] = 1.0  # m
    sys.constants['acc_due_to_grav'] = 9.81  # m/s**2

    sys.coordinates['angle'] = np.deg2rad(15.0)  # rad
    sys.speeds['angle_vel'] = 0.0  # rad/s

    assert list(sys.states.keys()) == ['angle', 'angle_vel']
    assert list(sys.states.values()) == [np.deg2rad(15.0), 0.0]

    def equations_of_motion(angle, angle_vel, length, acc_due_to_grav):

        l = length
        g = acc_due_to_grav

        thetad = angle_vel
        omegad = -g / l * np.sin(angle)

        return thetad, omegad

    # function that computes the right hand side of the equations of motion in
    # first order form
    sys.diff_eq_func = equations_of_motion

    def height(angle, length):
        """The Y coordinate points in the opposite of gravity, i.e. up. The X
        coordinate points to the right."""
        return -length * np.cos(angle)

    sys.add_measurement('height', height)

    def sway(angle, length):
        """The X coordinate points to the right."""
        return length * np.sin(angle)

    sys.add_measurement('sway', sway)

    def potential_energy(length, height, mass, acc_due_to_grav):
        return mass * acc_due_to_grav * (length + height)

    sys.add_measurement('potential_energy', potential_energy)

    def kinetic_energy(mass, length, angle_vel):
        return mass * (length * angle_vel)**2 / 2.0

    sys.add_measurement('kinetic_energy', kinetic_energy)

    def total_energy(kinetic_energy, potential_energy):
        return kinetic_energy + potential_energy

    sys.add_measurement('total_energy', total_energy)

    sample_rate = 30
    sys.free_response(5.0, sample_rate=sample_rate)

    def plot_config(mass, length, sway, height, time, sway__hist, height__hist,
                    time__hist, potential_energy__hist):
        fig, axes = plt.subplots(1, 2)

        circle_radius = mass / 2.0 / 10.0

        axes[0].set_xlim((-length - circle_radius, length + circle_radius))
        axes[0].set_ylim((-length - circle_radius, length + circle_radius))
        axes[0].set_title('Pendulum')
        axes[0].set_aspect('equal')
        xlabel = axes[0].set_xlabel('Time: {:.3f}'.format(time))

        path_lines = axes[0].plot(sway__hist, height__hist, color='red')
        rod_lines = axes[0].plot([0, sway], [0, height])

        circle = Circle((sway, height), radius=circle_radius)
        axes[0].add_patch(circle)

        axes[1].set_ylim((0, 0.4))
        axes[1].set_xlim((0, 5))
        axes[1].set_xlabel('Time [s]')
        axes[1].set_ylabel('Potential Energy [J]')
        pe_lines = axes[1].plot(time__hist, potential_energy__hist)

        plt.tight_layout()

        return fig, circle, rod_lines, path_lines, pe_lines, xlabel

    sys.config_plot_func = plot_config

    def update_plot(sway, height, time, time__hist, sway__hist, height__hist,
                    potential_energy__hist, circle, rod_lines, path_lines,
                    pe_lines, xlabel):
        circle.center = sway, height
        rod_lines[0].set_data([0, sway], [0, height])
        path_lines[0].set_data(sway__hist, height__hist)
        pe_lines[0].set_data(time__hist, potential_energy__hist)
        xlabel.set_text('Time: {:.3f}'.format(time))

        return circle, rod_lines, path_lines, pe_lines, xlabel

    sys.config_plot_update_func = update_plot


def test_sdof_trifilar_pendulum():

    sys = SingleDoFNonLinearSystem()

    sys.constants['m'] = 1  # kg
    sys.constants['r'] = 0.3  # m
    sys.constants['l'] = 0.75  # m
    sys.constants['g'] = 9.81  # m/s**2
    sys.constants['I'] = 0.3**2  # kg m**2

    sys.coordinates['theta'] = 0.2  # rad
    with pytest.raises(ValueError):
        sys.coordinates['omega'] = 0.0  # rad/s
    sys.speeds['omega'] = 0.0  # rad/s

    def eval_rhs(theta, omega, I, m, r, l, g):
        return theta

    with pytest.raises(ValueError):  # wrong num of return args in diff_eq_func
        sys.diff_eq_func = eval_rhs

    def eval_rhs(theta, omega, I, m, r, l, g):
        theta_dot = omega
        omega_dot = (-m*r**2*(2*g*(l**2 + 2*r**2*np.cos(theta) -
                                   2*r**2)**4*np.sin(theta) +
                     2*r**4*(l**2 + 2*r**2*np.cos(theta) -
                             2*r**2)**(5/2)*omega**2*np.sin(theta)**3 +
                     r**2*(l**2 + 2*r**2*np.cos(theta) -
                           2*r**2)**(7/2)*omega**2*np.sin(2*theta)) /
                     (2*(I*(l**2 + 2*r**2*np.cos(theta) - 2*r**2) +
                         m*r**4*np.sin(theta)**2)*(l**2 + 2*r**2*np.cos(theta) -
                                                   2*r**2)**(7/2)))
        return theta_dot, omega_dot

    sys.diff_eq_func = eval_rhs

    # should work with shape(m, 2n)
    res = sys._ode_eval_func(np.array([[0.1, 0.2],
                                       [0.1, 0.2],
                                       [0.1, 0.2]]), [0.1, 0.2, 0.3])
    assert res.shape == (3, 2)

    # reset these
    sys.coordinates['theta'] = 0.2  # rad
    sys.speeds['omega'] = 0.0  # rad/s
    sys._time['t'] = 0
    # should work with shape(m, 2n, 1)
    res = sys._ode_eval_func(np.random.random(3 * 2 * 1).reshape(3, 2, 1),
                             [0.1, 0.2, 0.3])
    assert res.shape == (3, 2, 1)

    # resest these
    sys.coordinates['theta'] = 0.2  # rad
    sys.speeds['omega'] = 0.0  # rad/s
    sys._time['t'] = 0
    traj = sys.free_response(2.0)

    # NOTE : See https://github.com/moorepants/resonance/issues/128
    sys.constants['l'] = 1.0  # m
    traj2 = sys.free_response(2.0)
    assert not traj.equals(traj2)

    assert 'theta' in traj.columns
    assert 'omega' in traj.columns
    assert 'theta_acc' in traj.columns

    desc = sys.__str__()
    expected_desc = """\
System name: SingleDoFNonLinearSystem

Differential equations function defined: True
Configuration plot function defined: False
Configuration update function defined: False

Constants
=========
m = 1.00000
r = 0.30000
l = 1.00000
g = 9.81000
I = 0.09000

Coordinates
===========
theta = 0.20000

Speeds
======
omega = d(theta)/dt = 0.00000

Measurements
============
"""
    assert desc == expected_desc


def test_multi_dof_nonlinear_system():

    sys = MultiDoFNonLinearSystem()

    sys.constants['m'] = 1.0  # kg
    sys.constants['k'] = 0.5  # N/m

    # NOTE : This could confuse the students because we couldn't do this with
    # linear systems. You had to pass these values into the response methods.
    sys.constants['Fo'] = 1.0  # N
    sys.constants['w'] = 3.0  # rad/s

    # NOTE : The order of declaration will define the order of the states.
    sys.coordinates['x2'] = 0.2  # m
    sys.coordinates['x1'] = 0.1  # m

    # should be in order of entry
    assert list(sys.coordinates.keys()) == ['x2', 'x1']

    # TODO : How do we know which speed corresponds to which coordinate
    # derivative?
    sys.speeds['v1'] = 0.01  # m/s
    sys.speeds['v2'] = 0.02  # m/s

    # should be in order of entry
    assert list(sys.speeds.keys()) == ['v1', 'v2']

    assert list(sys.states.keys()) == ['x2', 'x1', 'v1', 'v2']

    def rhs(x1, x2, v1, v2, m, k, Fo, w, time):
        # two masses connected by springs in series sliding on frictionless
        # surface with one spring attached to wall with sinusoidal forcing on
        # the end spring
        x1d = x1 * time
        x2d = x2 * time
        v1d = v1 * time
        v2d = v2 * time
        # NOTE : Order of output must match sys.states!
        return x2d, x1d, v1d, v2d

    sys.diff_eq_func = rhs

    desc = sys.__str__()
    expected_desc = """\
System name: MultiDoFNonLinearSystem

Differential equations function defined: True
Configuration plot function defined: False
Configuration update function defined: False

Constants
=========
m = 1.00000
k = 0.50000
Fo = 1.00000
w = 3.00000

Coordinates
===========
x2 = 0.20000
x1 = 0.10000

Speeds
======
v1 = d(x2)/dt = 0.01000
v2 = d(x1)/dt = 0.02000

Measurements
============
"""
    assert desc == expected_desc

    # should work with shape(2n,)
    x = np.random.random(4)
    t = 0.1
    res = sys._ode_eval_func(x, t)
    assert res.shape == (4,)
    np.testing.assert_allclose(res, x * t)
    sys.coordinates['x2'] = 0.2
    sys.coordinates['x1'] = 0.1
    sys.speeds['v1'] = 0.01
    sys.speeds['v2'] = 0.02
    sys._time['t'] = 0.0

    # should work with shape(1, 2n)
    x = np.random.random(4).reshape(1, 4)
    t = 0.1
    res = sys._ode_eval_func(x, t)
    assert res.shape == (1, 4)
    np.testing.assert_allclose(res, x * t)
    sys.coordinates['x2'] = 0.2
    sys.coordinates['x1'] = 0.1
    sys.speeds['v1'] = 0.01
    sys.speeds['v2'] = 0.02
    sys._time['t'] = 0.0

    # should work with shape(1, 2n, 1)
    x = np.random.random(4).reshape(1, 4, 1)
    t = 0.1
    res = sys._ode_eval_func(x, t)
    assert res.shape == (1, 4, 1)
    np.testing.assert_allclose(res, x * t)
    sys.coordinates['x2'] = 0.2
    sys.coordinates['x1'] = 0.1
    sys.speeds['v1'] = 0.01
    sys.speeds['v2'] = 0.02
    sys._time['t'] = 0.0

    # should work with shape(m, 2n)
    x = np.random.random(10 * 4).reshape(10, 4)
    t = np.random.random(10)
    res = sys._ode_eval_func(x, t)
    assert res.shape == (10, 4)
    np.testing.assert_allclose(res, x * t[:, np.newaxis])
    sys.coordinates['x2'] = 0.2
    sys.coordinates['x1'] = 0.1
    sys.speeds['v1'] = 0.01
    sys.speeds['v2'] = 0.02
    sys._time['t'] = 0.0

    # should work with shape(m, 2n, 1)
    x = np.random.random(10 * 4).reshape(10, 4, 1)
    t = np.random.random(10)
    res = sys._ode_eval_func(x, t)
    assert res.shape == (10, 4, 1)
    np.testing.assert_allclose(res, x * t[:, np.newaxis, np.newaxis])
    sys.coordinates['x2'] = 0.2
    sys.coordinates['x1'] = 0.1
    sys.speeds['v1'] = 0.01
    sys.speeds['v2'] = 0.02
    sys._time['t'] = 0.0

    # NOTE : Order of args does not matter here in the function signature.
    def rhs(x1, x2, v1, v2, m, k, Fo, w, time):
        # two masses connected by springs in series sliding on frictionless
        # surface with one spring attached to wall with sinusoidal forcing on
        # the end spring
        x1d = v1
        x2d = v2
        v1d = (-k * x1 + k * (x2 - x1)) / m
        v2d = (-k * (x2 - x1) + Fo * np.cos(w * time)) / m
        # NOTE : Order of output must match sys.states!
        return x2d, x1d, v1d, v2d

    sys.diff_eq_func = rhs

    traj = sys.free_response(5.0)

    for s in sys.states.keys():
        assert s in traj.columns

    # NOTE : The two integrators do not give the same answer, but they can be
    # compared to a low precision.
    traj = sys.free_response(3.0)
    traj2 = sys.free_response(3.0, integrator='lsoda', rtol=1.49012e-8)
    # TODO : This check it not working even at low precision. The results are
    # pretty different.
    #assert_frame_equal(traj, traj2, check_less_precise=1)


def test_measurements_in_diff_eq_func():
    # NOTE : Tests issue 127
    # https://github.com/moorepants/resonance/issues/127

    sys1 = SingleDoFNonLinearSystem()

    sys1.constants['m'] = 1
    sys1.constants['k'] = 100
    sys1.coordinates['x'] = 1
    sys1.speeds['v'] = 0

    def f(x, k):
        return k * x

    sys1.add_measurement('spring_force', f)

    def rhs(x, v, m, spring_force):
        xdot = v
        vdot = -spring_force / m
        return xdot, vdot

    sys1.diff_eq_func = rhs

    traj1 = sys1.free_response(2)
    del traj1['spring_force']

    sys2 = SingleDoFNonLinearSystem()

    sys2.constants['m'] = 1
    sys2.constants['k'] = 100
    sys2.coordinates['x'] = 1
    sys2.speeds['v'] = 0

    def rhs(x, v, m, k):
        xdot = v
        vdot = -k * x / m
        return xdot, vdot

    sys2.diff_eq_func = rhs

    traj2 = sys2.free_response(2)

    assert_frame_equal(traj1, traj2)
