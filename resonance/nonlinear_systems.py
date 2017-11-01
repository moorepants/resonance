from inspect import getargspec

import numpy as np
import scipy as sp
import scipy.integrate  # scipy doesn't import automatically
import matplotlib as mp

from .system import System as _System


class SingleDoFNonLinearSystem(_System):
    """This is the abstract base class for any single degree of freedom
    nonlinear system.  It can be sub-classed to make a custom system or the
    necessary methods can be added dynamically."""

    def __init__(self):

        super(SingleDoFNonLinearSystem, self).__init__()

        self._diff_eq_func = None

    @property
    def diff_eq_func(self):
        """A function that returns the time derivatives of the coordinate and
        speed, i.e. computes the right hand side of the two order differential
        equations. This equation looks like the following for linear motion:

            dy
            -- = f(t, y, p1, p2, ..., pO)
            dt

        where:

            - t: a time value
            - y: a vector containing the coordinate value and the speed value
            - p: any number of constants or measurements, O is the number of
              constants.

        Notes
        =====
        Your function should be able to operate on 1d arrays as inputs, i.e.
        use numpy math functions in your function, e.g. ``numpy.sin`` instead
        of ``math.sin``.

        Example
        =======
        >>> sys = SingleDoFNonLinearSystem()
        >>> sys.constants['gravity'] = 9.8  # m/s**2
        >>> sys.constants['length'] = 1.0  # m
        >>> sys.constnats['mass'] = 0.5  # kg
        >>> sys.coordinates['theta'] = 0.3  # rad
        >>> sys.speeds['omega'] = 0.0  # rad/s
        >>> def rhs(theta, omega, gravity, length, mass):
        >>>     # Represents a linear model of a simple pendulum:
        ...     #  m * l**2 ω' + m * g * l * sin(θ) = 0
        ...     thetad = omega
        ...     omegad = (-m*g*l*np.sin(theta)) / m / l**2
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

    @property
    def _array_rhs_eval_func(self):

        diff_eq_func_arg_names = getargspec(self.diff_eq_func).args
        diff_eq_func_args = [self._get_par_vals(k) for k in diff_eq_func_arg_names]

        coord_name = list(self.coordinates.keys())[0]
        speed_name = list(self.speeds.keys())[0]

        coord_idx = diff_eq_func_arg_names.index(coord_name)
        speed_idx = diff_eq_func_arg_names.index(speed_name)

        def eval_rhs(x, t):
            diff_eq_func_args[coord_idx] = x[0]
            diff_eq_func_args[speed_idx] = x[1]
            return np.asarray(self.diff_eq_func(*diff_eq_func_args))

        return eval_rhs

    def _integrate_equations_of_motion(self, times, integrator='rungakutta4'):

        x0 = list(self.coordinates.values())[0]
        v0 = list(self.speeds.values())[0]

        initial_conditions = np.array([x0, v0])

        method_name = '_integrate_with_{}'.format(integrator)
        integrator_method = getattr(self, method_name)

        return integrator_method(initial_conditions, times)

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
        return sp.integrate.odeint(self._array_rhs_eval_func,
                                   initial_conditions, times)

    def _integrate_with_rungakutta4(self, initial_conditions, times):
        """4th-order Runge-Kutta integration.

        Parameters
        ----------
        initial_conditions : array_like, shape(n,)
            Initial values of the state variables.
        times : array, shape(m,)
            Array of time values at which to solve.

        Returns
        -------
        x : ndarray, shape(m, n)
            Array containing the values of the state variables at the
            specified time values in `t`.

        """

        def _rk4(t, dt, x, f, args=None):
            """4th-order Runge-Kutta integration step."""
            x = np.asarray(x)
            if args is None:
                args = []
            k1 = np.asarray(f(x, t, *args))
            k2 = np.asarray(f(x + 0.5*dt*k1, t + 0.5*dt, *args))
            k3 = np.asarray(f(x + 0.5*dt*k2, t + 0.5*dt, *args))
            k4 = np.asarray(f(x + dt*k3, t + dt, *args))
            return x + dt*(k1 + 2*k2 + 2*k3 + k4)/6.0

        x = np.zeros((len(times), len(initial_conditions)))
        x[0, :] = initial_conditions
        for i in range(1, len(times)):
            dt = times[i] - times[i-1]
            x[i] = _rk4(times[i], dt, x[i-1], self._array_rhs_eval_func)
        return x

    def _generate_state_trajectories(self, times):
        """This method should return arrays for position, velocity, and
        acceleration of the coordinates."""
        int_res = self._integrate_equations_of_motion(times)

        res = self._array_rhs_eval_func(int_res.T, times)

        return int_res[:, 0], int_res[:, 1], res[1, :]


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
            The torsional Coulomb coefficient of friction at the pivot joint.
            The joint is clamped via a 1 Nm force.
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
                coeff_of_friction, acc_due_to_gravity):

            Irod_O = rod_mass * rod_length**2 / 3
            Ibob_P = bob_mass * bob_radius**2 / 2
            Ibob_O = Ibob_P + bob_mass * rod_length**2
            I = Irod_O + Ibob_O

            angle_dot = angle_vel
            angle_vel_dot = -(coeff_of_friction * np.sign(angle_vel) +
                              acc_due_to_gravity * rod_length * (bob_mass +
                              rod_mass / 2.0) * np.sin(angle)) / I

            # NOTE : These have to be in the correct order that matches
            # System.states, otherwise there is not way to detect which order
            # the user selected.
            return angle_dot, angle_vel_dot

            # TODO : Maybe we can let them use a dictionary as output to
            # specifically label things?
            #return {'angle': angle_dot, 'angle_vel': angle_vel_dot}

        self.diff_eq_func = rhs
