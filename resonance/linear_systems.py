import math
from inspect import getargspec
import warnings
import matplotlib.pyplot as plt
from matplotlib.patches import Circle

import numpy as np

from .system import System as _System


class SingleDoFLinearSystem(_System):

    @staticmethod
    def _natural_frequency(mass, stiffness):
        """Returns the real or complex valued natural frequency of the
        system."""
        wn = np.lib.scimath.sqrt(stiffness / mass)
        if isinstance(wn, complex):
            msg = ('The combination of system constants produces a complex '
                   'natural frequency, which results in an unstable system.')
            warnings.warn(msg)
        return wn

    @staticmethod
    def _damping_ratio(mass, damping, natural_frequency):
        zeta = damping / 2.0 / mass / natural_frequency
        if zeta * natural_frequency < 0.0:
            msg = ('The combination of system constants produces a negative '
                   'damping ratio, which results in an unstable system.')
            warnings.warn(msg)
        return zeta

    @staticmethod
    def _damped_natural_frequency(natural_frequency, damping_ratio):
        return natural_frequency * np.sqrt(1.0 - damping_ratio**2)

    def _solution_func(self):

        m, c, k = self._canonical_coefficients()
        omega_n = self._natural_frequency(m, k)

        if math.isclose(c, 0.0):
            if isinstance(omega_n, complex):
                sol_func = self._no_damping_unstable_solution
            else:
                sol_func = self._no_damping_solution
        else:  # damping, so check zeta
            zeta = self._damping_ratio(m, c, omega_n)
            if zeta < 1.0:
                sol_func = self._underdamped_solution
            elif zeta > 1.0:
                sol_func = self._overdamped_solution
            elif math.isclose(zeta, 1.0):
                sol_func = self._critically_damped_solution
            else:
                msg = 'No valid simulation solution with these parameters.'
                raise ValueError(msg)

        return sol_func

    def _no_damping_solution(self, time):

        t = time

        m, c, k = self._canonical_coefficients()

        wn = self._natural_frequency(m, k)

        x0 = list(self.coordinates.values())[0]
        v0 = list(self.speeds.values())[0]

        c1 = v0 / wn
        c2 = x0

        pos = c1 * np.sin(wn * t) + c2 * np.cos(wn * t)
        vel = c1 * wn * np.cos(wn * t) - c2 * wn * np.sin(wn * t)
        acc = -c1 * wn**2 * np.sin(wn * t) - c2 * wn**2 * np.cos(wn * t)

        return pos, vel, acc

    def _no_damping_unstable_solution(self, time):

        t = time

        m, c, k = self._canonical_coefficients()

        wn = self._natural_frequency(m, k).imag

        x0 = list(self.coordinates.values())[0]
        v0 = list(self.speeds.values())[0]

        # TODO : Verify these are correct.
        c1 = v0 / wn
        c2 = x0

        pos = c1 * np.sinh(wn * t) + c2 * np.cosh(wn * t)
        vel = wn * (c1 * np.cosh(wn * t) + c2 * np.sinh(wn * t))
        acc = wn**2 * (c1 * np.sinh(wn * t) + c2 * np.cosh(wn * t))

        return pos, vel, acc

    def _underdamped_solution(self, time):

        t = time

        m, c, k = self._canonical_coefficients()

        wn = self._natural_frequency(m, k)
        z = self._damping_ratio(m, c, wn)
        wd = self._damped_natural_frequency(wn, z)

        x0 = list(self.coordinates.values())[0]
        v0 = list(self.speeds.values())[0]

        A = np.sqrt(((v0 + z * wn * x0)**2 + (x0 * wd)**2) / wd**2)
        phi = np.arctan2(x0 * wd, v0 + z * wn * x0)

        pos = A * np.exp(-z * wn * t) * np.sin(wd * t + phi)

        vel = (A * -z * wn * np.exp(-z * wn * t) * np.sin(wd * t + phi) +
               A * np.exp(-z * wn * t) * wd * np.cos(wd * t + phi))

        acc = (A * (-z * wn)**2 * np.exp(-z * wn * t) * np.sin(wd * t + phi) +
               A * -z * wn * np.exp(-z * wn * t) * wd * np.cos(wd * t + phi) +
               A * -z * wn * np.exp(-z * wn * t) * wd * np.cos(wd * t + phi) -
               A * np.exp(-z * wn * t) * wd**2 * np.sin(wd * t + phi))

        return pos, vel, acc

    def _overdamped_solution(self, time):

        t = time

        m, c, k = self._canonical_coefficients()

        wn = self._natural_frequency(m, k)
        z = self._damping_ratio(m, c, wn)

        x0 = list(self.coordinates.values())[0]
        v0 = list(self.speeds.values())[0]

        a1 = ((-v0 + (-z + np.sqrt(z**2 - 1)) * wn * x0) / 2 / wn /
              np.sqrt(z**2 - 1))
        a2 = ((v0 + (z + np.sqrt(z**2 - 1)) * wn * x0) / 2 / wn /
              np.sqrt(z**2 - 1))

        time_const = wn * np.sqrt(z**2 - 1)

        pos = np.exp(-z*wn*t)*(a1*np.exp(-time_const*t) +
                               a2*np.exp(time_const*t))

        vel = (-z*wn*np.exp(-z*wn*t)*(a1*np.exp(-time_const*t) +
                                      a2*np.exp(time_const*t)) +
               np.exp(-z*wn*t)*(-a1*time_const*np.exp(-time_const*t) +
                                a2*time_const*np.exp(time_const*t)))

        acc = ((-z*wn)**2*np.exp(-z*wn*t)*(a1*np.exp(-time_const*t) +
                                           a2*np.exp(time_const*t)) +
               -z*wn*np.exp(-z*wn*t)*(-a1*time_const*np.exp(-time_const*t) +
                                      a2*time_const*np.exp(time_const*t)) +
               -z*wn*np.exp(-z*wn*t)*(-a1*time_const*np.exp(-time_const*t) +
                                      a2*time_const*np.exp(time_const*t)) +
               np.exp(-z*wn*t)*(a1*time_const**2*np.exp(-time_const*t) +
                                a2*time_const**2*np.exp(time_const*t)))

        return pos, vel, acc

    def _critically_damped_solution(self, time):

        t = time

        m, c, k = self._canonical_coefficients()

        wn = self._natural_frequency(m, k)

        x0 = list(self.coordinates.values())[0]
        v0 = list(self.speeds.values())[0]

        a1 = x0
        a2 = v0 + wn * x0

        pos = (a1 + a2 * t) * np.exp(-wn * t)
        vel = a2 * np.exp(-wn * t) + (a1 + a2 * t) * -wn * np.exp(-wn * t)
        acc = (a2 * -wn * np.exp(-wn * t) + a2 * -wn * np.exp(-wn * t) +
               (a1 + a2 * t) * wn**2 * np.exp(-wn * t))

        return pos, vel, acc

    def period(self):
        """Returns the (damped) period of oscillation of the coordinate in
        seconds."""
        m, c, k = self._canonical_coefficients()
        wn = self._natural_frequency(m, k)
        z = self._damping_ratio(m, c, wn)
        return 2.0 * np.pi / self._damped_natural_frequency(wn, z)

    def _generate_state_trajectories(self, times):

        sol_func = self._solution_func()

        return sol_func(times)


class BookOnCupSystem(SingleDoFLinearSystem):
    """This system represents dynamics of a typical engineering textbook set
    atop a cylinder (a coffee cup) such that the book can vibrate without slip
    on the curvature of the cup. It is described by:

    Attributes
    ==========
    constants
        thickness, t [meters]
            the thickness of the book
        length, l [meters]
            the length of the edge of the book which is tagent to the cup's
            surface
        mass, m [kilograms]
            the mass of the book
        radius, r [meters]
            the outer radius of the cup
    coordinates
        book_angle, theta [radians]
            the angle of the book with respect to the gravity vector
    speeds
        book_angle_vel, theta [radians]
            the angular rate of the book with repsect to the gravity vector

    """

    def __init__(self):

        super(BookOnCupSystem, self).__init__()

        self.constants['thickness'] = 0.029  # m
        self.constants['length'] = 0.238  # m
        self.constants['radius'] = 0.042  # m
        self.constants['mass'] = 1.058  # kg
        self.coordinates['book_angle'] = 0.0  # rad
        self.speeds['book_angle_vel'] = 0.0  # rad/s

    # TODO : This needs to be added in the super class with the add_coef_func()
    # method.
    def _canonical_coefficients(self):
        """A 1 DoF second order system should return the mass, damping, and
        stiffness coefficients."""

        def coeffs(thickness, length, radius):
            """Students will write this function themselves and pass it into
            the class via self.add_coef_func() when they get to modeling."""
            g = 9.81
            m = thickness**2 / 3 + length**2 / 12
            c = 0.0
            k = g * radius - g * thickness / 2
            return m, c, k

        args = [self._get_par_vals(k) for k in getargspec(coeffs).args]

        return coeffs(*args)


class TorsionalPendulumSystem(SingleDoFLinearSystem):
    """This system represents dynamics of a simple torsional pendulum in which
    the torsionally elastic member's axis is aligned with gravity and the axis
    of the torsion member passes through the mass center of an object attached
    to it's lower end. The top of the torsion rod is rigidly attached to the
    "ceiling". It is described by:

    Attributes
    ==========
    constants
        rotational_inertia, I [kg m**2]
            The moment of inertia of the object attached to the pendulum.
        torsional_damping, C [N s / m]
            The viscous linear damping coefficient which represents any energy
            disipation from things like air resistance, slip, etc.
        torsional_stiffness, K [N / m]
            The linear elastic stiffness coefficient of the torsion member,
            typically a round slender rod.
    coordinates
        torsional_angle, theta [rad]
    speeds
        torsional_angle_vel, theta_dot [rad / s]

    """

    def __init__(self):

        super(TorsionalPendulumSystem, self).__init__()

        self.constants['rotational_inertia'] = 0.0  # kg m^2
        self.constants['torsional_damping'] = 0.0  # Ns/m
        self.constants['torsional_stiffness'] = 0.0  # N/m

        # TODO : When a coordinate is added the speed should be automatically
        # added.
        self.coordinates['torsion_angle'] = 0.0
        self.speeds['torsion_angle_vel'] = 0.0

    def _canonical_coefficients(self):

        def coeffs(rotational_inertia, torsional_damping, torsional_stiffness):
            return rotational_inertia, torsional_damping, torsional_stiffness

        args = [self._get_par_vals(k) for k in getargspec(coeffs).args]

        return coeffs(*args)


class CompoundPendulumSystem(SingleDoFLinearSystem):
    """This system represents dynamics of a simple compound pendulum in which a
    rigid body is attached via a revolute joint to a fixed point. Gravity acts
    on the pendulum to bring it to an equilibrium state and there is no
    friction in the joint. It is described by:

    Attributes
    ==========
    constants
        pendulum_mass, m [kg]
            The mass of the compound pendulum.
        inertia_about_joint, i [kg m**2]
            The moment of inertia of the compound pendulum about the revolute
            joint.
        joint_to_mass_center, l [m]
            The distance from the revolute joint to the mass center of the
            compound pendulum.
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

        super(CompoundPendulumSystem, self).__init__()

        self.constants['pendulum_mass'] = 0.0  # kg
        self.constants['inertia_about_joint'] = 0.0  # kg m**2
        self.constants['joint_to_mass_center'] = 0.0  # m
        self.constants['acc_due_to_gravity'] = 0.0  # m / s**2

        # TODO : When a coordinate is added the speed should be automatically
        # added.
        self.coordinates['angle'] = 0.0
        self.speeds['angle_vel'] = 0.0

    def _canonical_coefficients(self):

        def coeffs(pendulum_mass, inertia_about_joint, joint_to_mass_center,
                   acc_due_to_gravity):
            m = pendulum_mass
            i = inertia_about_joint
            l = joint_to_mass_center
            g = acc_due_to_gravity

            return i, 0.0, m * g * l

        args = [self._get_par_vals(k) for k in getargspec(coeffs).args]

        return coeffs(*args)


class SimplePendulumSystem(SingleDoFLinearSystem):
    """This system represents dynamics of a simple pendulum in which a point
    mass is fixed on a massless pendulum arm of some length to a revolute
    joint. Gravity acts on the pendulum to bring it to an equilibrium state and
    there is no friction in the joint. It is described by:

    Attributes
    ==========
    constants
        pendulum_mass, m [kg]
            The mass of the compound pendulum.
        pendulum_length, l [m]
            The distance from the revolute joint to the point mass location.
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

        super(SimplePendulumSystem, self).__init__()

        self.constants['pendulum_mass'] = 0.0  # kg
        self.constants['pendulum_length'] = 0.0  # m
        self.constants['acc_due_to_gravity'] = 0.0  # m / s**2

        # TODO : When a coordinate is added the speed should be automatically
        # added.
        self.coordinates['angle'] = 0.0
        self.speeds['angle_vel'] = 0.0

    def _canonical_coefficients(self):

        def coeffs(pendulum_mass, pendulum_length, acc_due_to_gravity):
            m = pendulum_mass
            l = pendulum_length
            g = acc_due_to_gravity

            return m * l**2, 0.0, m * g * l

        args = [self._get_par_vals(k) for k in getargspec(coeffs).args]

        return coeffs(*args)


class ClockPendulumSystem(SingleDoFLinearSystem):
    """This system represents dynamics of a simple compound pendulum in which a
    rigid body is attached via a revolute joint to a fixed point. Gravity acts
    on the pendulum to bring it to an equilibrium state and there is no
    friction in the joint. It is described by:

    Attributes
    ==========
    constants
        pendulum_mass, m [kg]
            The mass of the compound pendulum.
        inertia_about_joint, i [kg m**2]
            The moment of inertia of the compound pendulum about the revolute
            joint.
        joint_to_mass_center, l [m]
            The distance from the revolute joint to the mass center of the
            compound pendulum.
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
        self.constants['viscous_damping'] = 0.0  # N s / m
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

        def plot_config(bob_radius, rod_length, bob_sway, bob_height, time):

            fig, ax = plt.subplots(1, 1)

            ax.set_xlim((-rod_length - bob_radius,
                         rod_length + bob_radius))
            ax.set_ylim((-rod_length - bob_radius, 0.0))
            ax.set_title('Pendulum')
            ax.set_aspect('equal')
            xlabel = ax.set_xlabel('Time: {:.2f}'.format(time))

            # NOTE : zorder ensures the patch is on top of the line.
            rod_lines = ax.plot([0, bob_sway], [0, bob_height], linewidth=6,
                                zorder=1)[0]

            circle = Circle((bob_sway, bob_height), radius=bob_radius,
                            color='red')
            circle.set_zorder(2)
            ax.add_patch(circle)

            return fig, circle, rod_lines, xlabel

        self.config_plot_func = plot_config

        def update_plot(bob_sway, bob_height, time, circle, rod_lines, xlabel):
            xlabel.set_text('Time: {:.2f}'.format(time))
            circle.center = bob_sway, bob_height
            rod_lines.set_data([0, bob_sway], [0, bob_height])

        self.config_plot_update_func = update_plot

    def _canonical_coefficients(self):

        def coeffs(bob_mass, bob_radius, rod_mass, rod_length, viscous_damping,
                   acc_due_to_gravity):

            Irod_O = rod_mass * rod_length**2 / 3
            Ibob_P = bob_mass * bob_radius**2 / 2
            Ibob_O = Ibob_P + bob_mass * rod_length**2

            I = Irod_O + Ibob_O
            C = viscous_damping * rod_length**2
            K = acc_due_to_gravity * rod_length * (bob_mass + rod_mass / 2.0)

            return I, C, K

        args = [self._get_par_vals(k) for k in getargspec(coeffs).args]

        return coeffs(*args)


class MassSpringDamperSystem(SingleDoFLinearSystem):
    """This system represents dynamics of a mass connected to a spring and
    damper (dashpot). The mass moves horizontally without friction and is acted
    on horizontally by the spring and damper in parallel. The system is
    described by:

    Attributes
    ==========
    constants
        mass, M [kg]
            The system mass.
        damping, C [kg / s]
            The viscous linear damping coefficient which represents any energy
            dissipation from things like air resistance, slip, etc.
        stiffness, K [N / m]
            The linear elastic stiffness of the spring.
    coordinates
        position, x [m]
    speeds
        velocity, x_dot [m / s]

    """

    def __init__(self):

        super(MassSpringDamperSystem, self).__init__()

        self.constants['mass'] = 1.0  # m
        self.constants['damping'] = 0.0  # kg/s
        self.constants['stiffness'] = 100  # N/m

        # TODO : When a coordinate is added the speed should be automatically
        # added.
        self.coordinates['position'] = 0.0
        self.speeds['velocity'] = 0.0

    def _canonical_coefficients(self):

        def coeffs(mass, damping, stiffness):
            return mass, damping, stiffness

        args = [self._get_par_vals(k) for k in getargspec(coeffs).args]

        return coeffs(*args)
