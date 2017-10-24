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

    def _integrate_equations_of_motion(self, times, integrator='rungakutta4'):

        x0 = list(self.coordinates.values())[0]
        v0 = list(self.speeds.values())[0]

        initial_conditions = np.array([x0, v0])

        method_name = '_integrate_with_{}'.format(integrator)
        integrator_method = getattr(self, method_name)

        return integrator_method(initial_conditions, times)

    def _integrate_with_lsoda(self, initial_conditions, times):
        """This method should return the integration results in the form of
        odeint."""

        rhs = self.equations_of_motion

        # NOTE: the first two args will always be the state and then time
        args = tuple([self._get_par_vals(k) for k in getargspec(rhs).args[2:]])

        return sp.integrate.odeint(rhs, initial_conditions, times, args=args)

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
            k1 = np.asarray(f(x, t, *args))
            k2 = np.asarray(f(x + 0.5*dt*k1, t + 0.5*dt, *args))
            k3 = np.asarray(f(x + 0.5*dt*k2, t + 0.5*dt, *args))
            k4 = np.asarray(f(x + dt*k3, t + dt, *args))
            return x + dt*(k1 + 2*k2 + 2*k3 + k4)/6.0

        func = self.equations_of_motion
        args = tuple([self._get_par_vals(k) for k in getargspec(func).args[2:]])
        x = np.zeros((len(times), len(initial_conditions)))
        x[0, :] = initial_conditions
        for i in range(1, len(times)):
            dt = times[i] - times[i-1]
            x[i] = _rk4(times[i], dt, x[i-1], func, args)
        return x

    def _generate_state_trajectories(self, times):
        """This method should return arrays for position, velocity, and
        acceleration of the coordinates."""
        int_res = self._integrate_equations_of_motion(times)
        f = self.equations_of_motion
        args = tuple([self._get_par_vals(k) for k in getargspec(f).args[2:]])
        res = np.asarray(f(int_res.T, times, *args))
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

        def rhs(state, times, bob_mass, bob_radius, rod_mass, rod_length,
                coeff_of_friction, acc_due_to_gravity):

            angle = state[0]
            angle_vel = state[1]

            Irod_O = rod_mass * rod_length**2 / 3
            Ibob_P = bob_mass * bob_radius**2 / 2
            Ibob_O = Ibob_P + bob_mass * rod_length**2
            I = Irod_O + Ibob_O

            angle_dot = angle_vel
            angle_vel_dot = -(coeff_of_friction * np.sign(angle_vel) +
                              acc_due_to_gravity * rod_length * (bob_mass +
                              rod_mass / 2.0) * np.sin(angle)) / I

            return [angle_dot, angle_vel_dot]

        self.equations_of_motion = rhs
