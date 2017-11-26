from inspect import getargspec

import numpy as np
import scipy as sp
import scipy.integrate  # scipy doesn't import automatically
import matplotlib as mp
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, Wedge

from .system import System as _System, _SingleDoFCoordinatesDict


class MultiDoFNonLinearSystem(_System):
    """This is the abstract base class for any single degree of freedom
    nonlinear system. It can be sub-classed to make a custom system or the
    necessary methods can be added dynamically."""

    def __init__(self):

        super(MultiDoFNonLinearSystem, self).__init__()

        self._diff_eq_func = None

    @property
    def diff_eq_func(self):
        """A function that returns the time derivatives of the coordinates and
        speeds, i.e. computes the right hand side of the explicit first order
        differential equations. This equation looks like the following for
        linear motion::

            dx
            -- = f(t, q1, ..., qn, u1, ..., un, p1, p2, ..., pO)
            dt

        where:

        - x: [q1, ..., qn, u1, ..., un], the "state vector"
        - t: a time value
        - q: the coordinates
        - u: the speeds
        - p: any number of constants, O is the number of constants

        Your function should be able to operate on 1d arrays as inputs, i.e.
        use numpy math functions in your function, e.g. ``numpy.sin`` instead
        of ``math.sin``. Besides the constants, coordinates, and speeds, there
        is a special variable ``time`` that you can use to give the current
        value of time inside your function.

        .. note:: The function has to return the derivatives of the states in
           the order of the ``state`` attribute.

        .. warning:: Do not use measurements as a function argument. This may
           cause causality issues and is not yet supported. You are unlikely to
           get a correct answer if you use a measurement in this function.

        Example
        =======
        >>> sys = SingleDoFNonLinearSystem()
        >>> sys.constants['gravity'] = 9.8  # m/s**2
        >>> sys.constants['length'] = 1.0  # m
        >>> sys.constants['mass'] = 0.5  # kg
        >>> sys.constants['omega_b'] = 0.1  # rad/s
        >>> sys.coordinates['theta'] = 0.3  # rad
        >>> sys.speeds['omega'] = 0.0  # rad/s
        >>> sys.states
        {'theta': 0.3, 'omega': 0.0}  # note the order!
        >>> def rhs(theta, omega, gravity, length, mass, omega_b, time):
        ...     # Represents a linear model of a simple pendulum under
        ...     # sinusoidal torquing.
        ...     #  m * l**2 ω' + m * g * l * sin(θ) = sin(ω_b * t)
        ...     thetad = omega
        ...     omegad = (np.sin(omega_b * time) -
        ...               m*g*l*np.sin(theta)) / m / l**2
        ...     return thetad, omegad  # in order of sys.states
        >>> sys.diff_eq_func = rhs

        """
        return self._diff_eq_func

    @diff_eq_func.setter
    def diff_eq_func(self, func):
        self._measurements._check_for_duplicate_keys()
        # NOTE : This will throw an error if the function's args are not in the
        # system.
        [self._get_par_vals(k) for k in getargspec(func).args]
        self._diff_eq_func = func
        self._check_diff_eq()
        self._ode_eval_func = self._generate_array_rhs_eval_func()

    def _check_diff_eq(self):

        arg_names = getargspec(self.diff_eq_func).args
        arg_vals = [self._get_par_vals(k) for k in arg_names]

        msg = ('Your diff_eq_func does not return the correct number of '
               'state derivatives. Make sure the number and order of the '
               'states match the derivatives of the states you return.')

        res = self.diff_eq_func(*arg_vals)
        try:
            len(res)
        except TypeError:  # returns a single value, must return at least 2
            raise ValueError(msg)
        else:
            if len(res) != len(self.states):
                raise ValueError(msg)

    def _generate_array_rhs_eval_func(self):

        arg_names = getargspec(self.diff_eq_func).args

        coord_names = list(self.coordinates.keys())
        speed_names = list(self.speeds.keys())

        n = len(speed_names)

        def eval_rhs(x, t):
            # x is either shape(1, 2n), shape(m, 2n), shape(1, 2n, 1) or shape(m, 2n, 1)
            # t is either float or shape(m,)
            # TODO : This could be slow for large # coords/speeds.
            for i, cname in enumerate(coord_names):
                self.coordinates[cname] = x[:, i, ...]
            for i, sname in enumerate(speed_names):
                self.speeds[sname] = x[:, i + n, ...]

            if len(x.shape) == 3 and x.shape[-1] == 1:
                self._time['t'] = np.atleast_2d(t).T
            else:
                self._time['t'] = np.asarray(t)

            # this calculates the measurements
            arg_vals = [self._get_par_vals(k) for k in arg_names]

            # TODO : Would be nice not to have to create this every eval.
            x_dot = np.zeros_like(x)
            for i, dot in enumerate(self.diff_eq_func(*arg_vals)):
                x_dot[:, i, ...] = dot

            return x_dot

        return eval_rhs

    def _integrate_equations_of_motion(self, times, integrator='rungakutta4'):
        # TODO : This overrides the integrator option. Remove this once the
        # other integrator(s) work.
        integrator = 'rungakutta4'

        x0 = list(self.coordinates.values())
        v0 = list(self.speeds.values())

        if len(x0) != len(v0):
            msg = ('There is not and equal number of coordinates and speeds. '
                   'Make sure you have added one speed for each coordinate.')
            raise ValueError(msg)

        initial_conditions = np.hstack((x0, v0))

        method_name = '_integrate_with_{}'.format(integrator)
        integrator_method = getattr(self, method_name)

        # make sure rhs is up-to-date
        self._ode_eval_func = self._generate_array_rhs_eval_func()

        try:
            traj = integrator_method(initial_conditions, times)
        except:
            raise
        else:  # integration succeeds
            return traj

    def _integrate_with_lsoda(self, initial_conditions, times):
        """This method should return the integration results in the form of
        odeint.

        Parameters
        ==========
        initial_conditions : ndarray, shape(n,)
            The initial condition of each state.
        times : ndarray, shape(m,)
            The monotonically increasing time values.

        """
        return sp.integrate.odeint(self._ode_eval_func,
                                   initial_conditions, times)

    def _integrate_with_rungakutta4(self, initial_conditions, times):
        """4th-order Runge-Kutta integration.

            Array of time values at which to solve.

        Returns
        -------
        x : ndarray, shape(m, n)
            Array containing the values of the state variables at the
            specified time values in `t`.

        """

        def _rk4(t, dt, x, f, args=None):
            """4th-order Runge-Kutta integration step."""
            # x can have shape(2n, 1)
            # f returns shape(2n, 1)
            # t is a float
            # dt is a float
            x = x[np.newaxis, ...]
            if args is None:
                args = []
            k1 = f(x, t, *args)
            k2 = f(x + 0.5*dt*k1, t + 0.5*dt, *args)
            k3 = f(x + 0.5*dt*k2, t + 0.5*dt, *args)
            k4 = f(x + dt*k3, t + dt, *args)
            return x + dt*(k1 + 2*k2 + 2*k3 + k4)/6.0

        # m x 2n x 1
        x = np.zeros((len(times), len(initial_conditions), 1))
        x[0, :, 0] = initial_conditions
        for i in range(1, len(times)):
            dt = times[i] - times[i-1]
            # x[i] is 2n x 1
            x[i] = _rk4(times[i], dt, x[i-1], self._ode_eval_func)
        return x

    def _generate_state_trajectories(self, times, integrator='rungakutta4'):
        """This method should return arrays for position, velocity, and
        acceleration of the coordinates."""

        # m : num time samples
        # n : num coordinates/speeds

        # store values before integration
        coords = self.coordinates.copy()
        speeds = self.speeds.copy()
        time = self._time['t']

        try:
            # rows correspond to time, columns to states (m x 2n x 1)
            int_res = self._integrate_equations_of_motion(
                times, integrator=integrator)

            if int_res.shape != (len(times), len(self.states), 1):
                msg = ('Shape of trajectory from integration does not have '
                       'the correct shape.')
                raise ValueError(msg)

            # calculate the accelerations
            res = self._ode_eval_func(int_res, times)

            if res.shape != (len(times), len(self.states), 1):
                msg = ('Shape of derivatives does not have the correct shape.')
                raise ValueError(msg)
        except:
            raise
        finally:  # make sure to reset coords, speeds, time if anything fails
            # reset to values before integration
            for k, v in coords.items():
                self.coordinates[k] = v
            for k, v in speeds.items():
                self.speeds[k] = v
            self._time['t'] = time

        num_coords = len(self.coordinates)
        num_speeds = len(self.speeds)

        if num_coords != num_speeds:
            msg = ('You do not have the same number of coordinates as you do '
                   'speeds. There should be one speed for each coordinate.')
            raise ValueError(msg)

        pos = int_res[:, :num_coords, 0].T  # n x m
        vel = int_res[:, num_coords:, 0].T  # n x m
        acc = res[:, num_coords:, 0].T  # n x m

        return pos, vel, acc


class SingleDoFNonLinearSystem(MultiDoFNonLinearSystem):

    def __init__(self):

        super(SingleDoFNonLinearSystem, self).__init__()

        self._coordinates = _SingleDoFCoordinatesDict({})
        self._speeds = _SingleDoFCoordinatesDict({})
        self._measurements._coordinates = self._coordinates
        self._measurements._speeds = self._speeds


class ClockPendulumSystem(SingleDoFNonLinearSystem):
    """This system represents dynamics of a compound pendulum representing a
    clock pendulum. It is made up of a thin long cylindrical rod with a thin
    disc bob on the end. Gravity acts on the pendulum to bring it to an
    equilibrium state and there is option Coulomb friction in the joint. It is
    described by:

    Attributes
    ==========
    constants
        bob_mass, m_b [kg]
            The mass of the bob (a thin disc) on the end of the pendulum.
        bob_radius, r [m]
            The radius of the bob (a thin disc) on the end of the pendulum.
        rod_mass, m_r [kg]
            The mass of the then cylindrical rod.
        rod_length, l [m]
            The length of the rod which connects the pivot joint to the center
            of the bob.
        coeff_of_friction, mu [unitless]
            The Coulomb coefficient of friction between the materials of the
            pivot joint.
        joint_friction_radius, R [m]
            The radius of the contact disc at the pivot joint. The joint is
            assumed to be two flat discs pressed together.
        joint_clamp_force, F_N [N]
            The clamping force pressing the two flat discs together at the
            pivot joint.
        acc_due_to_gravity, g [m/s**2]
            The acceleration due to gravity.
    coordinates
        angle, theta [rad]
            The angle of the pendulum relative to the direction of gravity.
            When theta is zero the pendulum is hanging down in it's equilibrium
            state.
    speeds
        angle_vel, theta_dot [rad / s]
            The angular velocity of the pendulum about the revolute joint axis.

    """

    def __init__(self):

        super(ClockPendulumSystem, self).__init__()

        self.constants['bob_mass'] = 0.1  # kg
        self.constants['bob_radius'] = 0.03  # m
        self.constants['rod_mass'] = 0.1  # kg
        self.constants['rod_length'] = 0.2799  # m
        self.constants['coeff_of_friction'] = 0.0  # unitless
        self.constants['joint_friction_radius'] = 0.03  # m
        self.constants['joint_clamp_force'] = 1.0  # N
        self.constants['acc_due_to_gravity'] = 9.81  # m / s**2

        self.coordinates['angle'] = 0.0
        self.speeds['angle_vel'] = 0.0

        def bob_height(angle, rod_length):
            """The Y coordinate of the bob. The Y coordinate points in the
            opposite of gravity, i.e. up. The X coordinate points to the
            right."""
            return -rod_length * np.cos(angle)

        self.add_measurement('bob_height', bob_height)

        def bob_sway(angle, rod_length):
            """The X coordinate of the bob center. The X coordinate points to
            the right."""
            return rod_length * np.sin(angle)

        self.add_measurement('bob_sway', bob_sway)

        def kinetic_energy(bob_mass, bob_radius, rod_length, bob_height,
                           rod_mass, angle_vel):

            v_bob = rod_length * angle_vel
            v_rod = rod_length / 2 * angle_vel

            I_bob = bob_mass * bob_radius**2 / 2
            I_rod = rod_mass * rod_length**2 / 12

            KE_bob = bob_mass * v_bob**2 / 2 + I_bob * angle_vel**2 / 2
            KE_rod = rod_mass * v_rod**2 / 2 + I_rod * angle_vel**2 / 2

            return KE_rod + KE_bob

        self.add_measurement('kinetic_energy', kinetic_energy)

        def potential_energy(bob_mass, rod_mass, rod_length, bob_height,
                             acc_due_to_gravity, angle):
            PE_bob = bob_mass * acc_due_to_gravity * (rod_length - rod_length *
                                                      np.cos(angle))
            PE_rod = rod_mass * acc_due_to_gravity * (rod_length / 2 -
                                                      rod_length / 2 *
                                                      np.cos(angle))
            return PE_bob + PE_rod

        self.add_measurement('potential_energy', potential_energy)

        def total_energy(kinetic_energy, potential_energy):
            return kinetic_energy + potential_energy

        self.add_measurement('total_energy', total_energy)

        def plot_config(bob_radius, rod_length, bob_sway, bob_height, time):

            fig, ax = mp.pyplot.subplots(1, 1)

            ax.set_xlim((-rod_length - bob_radius,
                         rod_length + bob_radius))
            ax.set_ylim((-rod_length - bob_radius, 0.0))
            ax.set_title('Pendulum')
            ax.set_aspect('equal')
            xlabel = ax.set_xlabel('Time: {:.2f}'.format(time))

            # NOTE : zorder ensures the patch is on top of the line.
            rod_lines = ax.plot([0, bob_sway], [0, bob_height], linewidth=6,
                                zorder=1)[0]

            circle = mp.patches.Circle((bob_sway, bob_height),
                                       radius=bob_radius, color='red')
            circle.set_zorder(2)
            ax.add_patch(circle)

            return fig, circle, rod_lines, xlabel

        self.config_plot_func = plot_config

        def update_plot(bob_sway, bob_height, time, circle, rod_lines, xlabel):
            xlabel.set_text('Time: {:.2f}'.format(time))
            circle.center = bob_sway, bob_height
            rod_lines.set_data([0, bob_sway], [0, bob_height])

        self.config_plot_update_func = update_plot

        def rhs(angle, angle_vel, bob_mass, bob_radius, rod_mass, rod_length,
                coeff_of_friction, joint_friction_radius, joint_clamp_force,
                acc_due_to_gravity):

            Irod_O = rod_mass * rod_length**2 / 3
            Ibob_P = bob_mass * bob_radius**2 / 2
            Ibob_O = Ibob_P + bob_mass * rod_length**2
            I = Irod_O + Ibob_O

            friction_torque = (2 / 3 * joint_friction_radius *
                               coeff_of_friction * joint_clamp_force *
                               np.sign(angle_vel))

            angle_dot = angle_vel
            angle_vel_dot = -(friction_torque +
                              acc_due_to_gravity * rod_length *
                              (bob_mass + rod_mass / 2.0) * np.sin(angle)) / I

            # NOTE : These have to be in the correct order that matches
            # System.states, otherwise there is no way to detect which order
            # the user selected.
            return angle_dot, angle_vel_dot

            # TODO : Maybe we can let them use a dictionary as output to
            # specifically label things?
            #return {'angle': angle_dot, 'angle_vel': angle_vel_dot}

        self.diff_eq_func = rhs


class BallChannelPendulumSystem(MultiDoFNonLinearSystem):

    def __init__(self):

        super(BallChannelPendulumSystem, self).__init__()

        self.constants['mp'] = 12/1000  # kg
        self.constants['mb'] = 3.5/1000  # kg
        self.constants['r'] = 0.1  # m
        self.constants['l'] = 0.2  # m
        self.constants['g'] = 9.81  # m/s**2

        self.coordinates['theta'] = np.deg2rad(10)
        self.coordinates['phi'] = np.deg2rad(-10)

        self.speeds['alpha'] = 0.0
        self.speeds['beta'] = 0.0

        def pend_y(l, theta):
            return (l - l * np.cos(theta))

        def pend_x(l, theta):
            return l * np.sin(theta)

        self.add_measurement('pend_x', pend_x)
        self.add_measurement('pend_y', pend_y)

        def ball_y(l, r, theta, phi):
            return l + r * np.cos(theta) - r * np.cos(theta + phi)

        def ball_x(l, r, theta, phi):
            return -r * np.sin(theta) + r * np.sin(theta + phi)

        self.add_measurement('ball_x', ball_x)
        self.add_measurement('ball_y', ball_y)

        def trough_x(r, theta):
            return -r * np.sin(theta)

        def trough_y(l, r, theta):
            return l + r * np.cos(theta)

        self.add_measurement('trough_x', trough_x)
        self.add_measurement('trough_y', trough_y)

        def create_plot(pend_x, pend_y, ball_x, ball_y,
                        trough_x, trough_y, l, r):
            # create a blank figure and set basic settings on the axis
            fig, ax = plt.subplots(1, 1)
            ax.set_xlim((-1, 1.0))
            ax.set_ylim((-r, l + r + r))
            ax.set_xlabel('x [m]')
            ax.set_ylabel('y [m]')
            ax.set_aspect('equal')

            ax.plot([0, 0], [0, l])

            pend_line = ax.plot([0, pend_x], [l, pend_y], color='red')[0]

            trough = Wedge((trough_x, trough_y), r, 180, 360, width=0.01)

            # circles are created by supplying an (x, y) pair and the radius
            ball = Circle((ball_x, ball_y), radius=0.02, color='black')
            bob = Circle((pend_x, pend_y), radius=0.05)

            ax.add_patch(trough)
            ax.add_patch(ball)
            ax.add_patch(bob)

            return fig, ball, bob, trough, pend_line

        self.config_plot_func = create_plot

        def update(pend_x, pend_y, ball_x, ball_y, l, theta, trough_x,
                   trough_y, ball, bob, trough, pend_line):
            ball.center = (ball_x, ball_y)
            bob.center = (pend_x, pend_y)
            pend_line.set_data([0, pend_x], [l, pend_y])
            trough.set_theta1(180 + np.rad2deg(theta))
            trough.set_theta2(360 + np.rad2deg(theta))
            trough.set_center((trough_x, trough_y))

        self.config_plot_update_func = update
