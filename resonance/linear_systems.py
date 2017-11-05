import math
from inspect import getargspec
import warnings
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, Rectangle

import numpy as np

from .system import System as _System
from .system import _SingleDoFCoordinatesDict


class SingleDoFLinearSystem(_System):
    """This is the abstract base class for any single degree of freedom linear
    system. It can be sub-classed to make a custom system or the necessary
    methods can be added dynamically."""

    def __init__(self):

        super(SingleDoFLinearSystem, self).__init__()

        self._coordinates = _SingleDoFCoordinatesDict({})
        self._speeds = _SingleDoFCoordinatesDict({})
        self._measurements._coordinates = self._coordinates
        self._measurements._speeds = self._speeds

        self._canonical_coeffs_func = None

    def _initial_conditions(self):
        x0 = list(self.coordinates.values())[0]
        v0 = list(self.speeds.values())[0]
        return x0, v0

    @property
    def canonical_coeffs_func(self):
        """A function that returns the three linear coefficients of the left
        hand side of a canonical second order ordinary differential equation.
        This equation looks like the following for linear motion:

            mv' + cv + kx = F(t)

        and like the following for angular motion:

            Iω' + cω + kθ = T(t)

        where:

            - m: mass of the moving particle
            - I: moment of inertia of a rigid body
            - c: viscous damping coefficient (linear or angular)
            - k: spring stiffness (linear or angular)
            - x: the positional coordinate of the mass
            - v: the positional speed of the mass
            - θ: the angular coordinate of the body
            - ω: the angular speed of the body

        The coefficients (m, c, k, I) must be defined in terms of the system's
        constants.

        Example
        =======
        >>> sys = SingleDoFLinearSystem()
        >>> sys.constants['gravity'] = 9.8  # m/s**2
        >>> sys.constants['length'] = 1.0  # m
        >>> sys.constnats['mass'] = 0.5  # kg
        >>> sys.coordinates['theta'] = 0.3  # rad
        >>> sys.speeds['omega'] = 0.0  # rad/s
        >>> def coeffs(gravity, length, mass):
        >>>     # Represents a linear model of a simple pendulum:
        ...     #  m * l**2 ω' + m * g * l * θ = 0
        ...     I = mass * length**2
        ...     c = 0.0
        ...     k = mass * gravity * length
        ...     return I, c, k
        >>> sys.canonical_coeffs_func = coeffs

        """
        return self._canonical_coeffs_func

    @canonical_coeffs_func.setter
    def canonical_coeffs_func(self, func):
        self._measurements._check_for_duplicate_keys()
        for k in getargspec(func).args:
            # NOTE : Measurements do not have to be time varying.
            if k not in (list(self.constants.keys()) +
                         list(self.measurements.keys())):
                msg = ('The function argument {} is not in constants or '
                       'measurements. Redefine your function in terms of '
                       'non-time varying parameters.')
                raise ValueError(msg.format(k))
        self._canonical_coeffs_func = func

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

    def _damped_natural_frequency(self, natural_frequency, damping_ratio):
        typ = self._solution_type()
        if typ == 'underdamped':
            return natural_frequency * np.sqrt(1.0 - damping_ratio**2)
        elif typ == 'no_damping_unstable' or typ == 'no_damping':
            return natural_frequency
        else:
            return 0.0

    def _normalized_form(self, m, c, k):
        wn = self._natural_frequency(m, k)
        z = self._damping_ratio(m, c, wn)
        wd = self._damped_natural_frequency(wn, z)
        return wn, z, wd

    def _solution_type(self):

        m, c, k = self._canonical_coefficients()
        omega_n = self._natural_frequency(m, k)

        if math.isclose(c, 0.0):
            if isinstance(omega_n, complex):
                return 'no_damping_unstable'
            else:
                return 'no_damping'
        else:  # damping, so check zeta
            zeta = self._damping_ratio(m, c, omega_n)
            if zeta < 1.0:
                return 'underdamped'
            elif math.isclose(zeta, 1.0):
                return 'critically_damped'
            elif zeta > 1.0:
                return 'overdamped'
            else:
                msg = 'No valid simulation solution with these parameters.'
                raise ValueError(msg)

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
            elif math.isclose(zeta, 1.0):
                sol_func = self._critically_damped_solution
            elif zeta > 1.0:
                sol_func = self._overdamped_solution
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

    def _damped_sinusoid(self, A, phi, t):
        """pos = A * exp(-z*wn*t) * sin(wd*t + phi)

        A and phi are different and depend on the particular solution.

        """

        m, c, k = self._canonical_coefficients()
        wn, z, wd = self._normalized_form(m, c, k)

        pos = A * np.exp(-z*wn*t) * np.sin(wd*t + phi)

        vel = (A * -z * wn * np.exp(-z*wn*t) * np.sin(wd*t + phi) +
               A * np.exp(-z*wn*t) * wd * np.cos(wd*t + phi))

        acc = (A * (-z * wn)**2 * np.exp(-z*wn*t) * np.sin(wd*t + phi) +
               A * -z * wn * np.exp(-z*wn*t) * wd * np.cos(wd*t + phi) +
               A * -z * wn * np.exp(-z*wn*t) * wd * np.cos(wd*t + phi) -
               A * np.exp(-z*wn*t) * wd**2 * np.sin(wd*t + phi))

        return pos, vel, acc

    def _underdamped_solution(self, time):

        t = time

        m, c, k = self._canonical_coefficients()
        wn, z, wd = self._normalized_form(m, c, k)

        x0, v0 = self._initial_conditions()

        A = np.sqrt(((v0 + z * wn * x0)**2 + (x0 * wd)**2) / wd**2)

        phi = np.arctan2(x0 * wd, v0 + z * wn * x0)

        pos, vel, acc = self._damped_sinusoid(A, phi, t)

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
                                           a2*np.exp(time_const*t)) -
               z*wn*np.exp(-z*wn*t)*(-a1*time_const*np.exp(-time_const*t) +
                                     a2*time_const*np.exp(time_const*t)) -
               z*wn*np.exp(-z*wn*t)*(-a1*time_const*np.exp(-time_const*t) +
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
        if z < 1.0:  # underdamped, no damping, or unstable
            return 2.0 * np.pi / self._damped_natural_frequency(wn, z)
        else:
            return np.inf

    def _generate_state_trajectories(self, times):

        sol_func = self._solution_func()

        return sol_func(times)

    def _canonical_coefficients(self):
        if self.canonical_coeffs_func is None:
            msg = ('There is no function available to calculate the canonical'
                   ' coeffcients.')
            raise ValueError(msg)
        else:
            f = self.canonical_coeffs_func
            args = [self._get_par_vals(k) for k in getargspec(f).args]
            return f(*args)

    def _periodic_forcing_steady_state(self, a0, an, bn, wT, t):

        M = t.shape[0]

        # scalars
        m, c, k = self._canonical_coefficients()
        wn, z, wd = self._normalized_form(m, c, k)

        N = an.shape[0]

        # column array of n values, shape(N, 1)
        n = np.arange(1, N+1)[:, np.newaxis]
        assert n.shape == (N, 1)

        # phase shift of each term in the series, shape(N, 1)
        theta_n = np.arctan2(2*z*wn*n*wT, wn**2-(n*wT)**2)
        assert theta_n.shape == (N, 1)

        # an is a col and t is 1D, so each row of xcn is a term
        # in the series at all times in t

        # shape(N, 1)
        denom = m * np.sqrt((wn**2 - (n*wT)**2)**2 + (2*z*wn*n*wT)**2)
        assert denom.shape == (N, 1)

        # shape(N, M)
        cwT = np.cos(n*wT*t - theta_n)
        swT = np.sin(n*wT*t - theta_n)
        assert cwT.shape == (N, M)
        assert swT.shape == (N, M)

        # shape(N, M)
        xcn = an / denom * cwT
        vcn = -an * n * wT / denom * swT
        acn = -an * (n * wT)**2 / denom * cwT
        assert xcn.shape == (N, M)

        # shape(N, M)
        xsn = bn / denom * swT
        vsn = bn * n * wT / denom * cwT
        asn = -bn * (n * wT)**2 / denom * swT
        assert xsn.shape == (N, M)

        # steady state solution (particular solution)
        # x is the sum of each xcn term (the rows)
        xss = a0 / 2 / k + np.sum(xcn + xsn, axis=0)
        vss = np.sum(vcn + vsn, axis=0)
        ass = np.sum(acn + asn, axis=0)
        assert xss.shape == (M, )

        return xss, vss, ass, n, theta_n, denom

    def _periodic_forcing_transient_A_phi(self, wT, n, a0, an, bn, theta_n,
                                          denom, t):
        # scalars
        m, c, k = self._canonical_coefficients()
        wn, z, wd = self._normalized_form(m, c, k)

        # the transient solution (homogeneous)
        x0, v0 = self._initial_conditions()

        c1 = np.sum((-np.sin(theta_n)*bn + np.cos(theta_n)*an) / denom)
        c2 = wT * np.sum((np.sin(theta_n)*an + np.cos(theta_n)*bn) * n / denom)

        phi = np.arctan2(wd*(2*c1*k+a0-2*k*x0),
                         2*c1*k*wn*z + 2*c2*k + a0*wn*z - 2*k*wn*x0*z -
                         2*k*v0)
        A = (-a0 / 2 + k * (-c1 + x0)) / k / np.sin(phi)

        return A, phi

    def periodic_forcing_response(self, twice_avg, cos_coeffs, sin_coeffs,
                                  frequency, final_time, initial_time=0.0,
                                  sample_rate=100,
                                  col_name='forcing_function'):
        """Returns the trajectory of the system's coordinates, speeds,
        accelerations, and measurements if a periodic forcing function defined
        by a Fourier series is applied as a force or torque in the same
        direction as the system's coordinate. The forcing function is defined
        as::

                                    N
            F(t) or T(t) = a0 / 2 + ∑ (an * cos(n*ω*t) + bn * sin(n*ω*t))
                                   n=1

        Where a0, a1...an, and b1...bn are the Fourier coefficients. If N=∞
        then the Fourier series can describe any periodic function with a
        period (2*π)/ω.

        Parameters
        ==========
        twice_avg : float
            Twice the average value over one cycle, a0.
        cos_coeffs : float or sequence of floats
            The N cosine Fourier coefficients: a1, ..., aN.
        sin_coeffs : float or sequence of floats
            The N sine Fourier coefficients: b1, ..., bN.
        frequency : float
            The frequency, ω, in radians per second corresponding to one full
            cycle of the function.
        final_time : float
            A value of time in seconds corresponding to the end of the
            simulation.
        initial_time : float, optional
            A value of time in seconds corresponding to the start of the
            simulation.
        sample_rate : integer, optional
            The sample rate of the simulation in Hertz (samples per second).
            The time values will be reported at the initial time and final
            time, i.e. inclusive, along with times space equally based on the
            sample rate.
        col_name : string, optional
            A valid Python identifier that will be used as the column name for
            the forcing function trajectory in the returned data frame.

        Returns
        =======
        pandas.DataFrame
            A data frame indexed by time with all of the coordinates and
            measurements as columns.

        """
        if self._solution_type() != 'underdamped':
            msg = 'Currently, only supported for underdamped systems.'
            raise ValueError(msg)

        # shape(N, 1)
        an = np.atleast_2d(cos_coeffs).T
        bn = np.atleast_2d(sin_coeffs).T

        if an.shape[1] != 1 or bn.shape[1] != 1:
            msg = 'an and bn must be 1D sequences or a single float.'
            raise ValueError(msg)

        if an.shape != bn.shape:
            raise ValueError('an and bn must be the same length')

        # shape (M,), M: number of time samples
        t = self._calc_times(final_time, initial_time, sample_rate)

        # scalars
        m, c, k = self._canonical_coefficients()
        wn, z, wd = self._normalized_form(m, c, k)
        wT = frequency

        xss, vss, ass, n, theta_n, denom = \
            self._periodic_forcing_steady_state(twice_avg, an, bn, wT, t)

        # the transient solution (homogeneous)

        A, phi = self._periodic_forcing_transient_A_phi(wT, n, twice_avg, an,
                                                        bn, theta_n, denom, t)
        xh, vh, ah = self._damped_sinusoid(A, phi, t)

        self.result = self._state_traj_to_dataframe(t, xh + xss, vh + vss,
                                                    ah + ass)

        if not col_name.isidentifier():
            msg = "'{}' is not a valid Python identifier."
            raise ValueError(msg.format(col_name))
        elif col_name in self.result.columns:
            raise ValueError('{} already taken.'.format(col_name))
        else:
            self.result[col_name] = twice_avg / 2 + \
                np.sum(an * n * np.cos(frequency * n * t) +
                       bn * n * np.sin(frequency * n * t), axis=0)

        return self.result

    def sinusoidal_forcing_response(self, amplitude, frequency, final_time,
                                    initial_time=0.0, sample_rate=100,
                                    col_name='forcing_function'):
        """Returns the trajectory of the system's coordinates, speeds,
        accelerations, and measurements if a sinusoidal forcing (or torquing)
        function defined by:

        F(t) = Fo * cos(ω * t)

        or

        T(t) = To * cos(ω * t)

        is applied to the moving body in the direction of the system's
        coordinate.

        Parameters
        ==========
        amplitude : float
            The amplitude of the forcing/torquing function, Fo or To, in
            Newtons or Newton-Meters.
        frequency : float
            The frequency, ω, in radians per second of the sinusoidal forcing.
        final_time : float
            A value of time in seconds corresponding to the end of the
            simulation.
        initial_time : float, optional
            A value of time in seconds corresponding to the start of the
            simulation.
        sample_rate : integer, optional
            The sample rate of the simulation in Hertz (samples per second).
            The time values will be reported at the initial time and final
            time, i.e. inclusive, along with times space equally based on the
            sample rate.
        col_name : string, optional
            A valid Python identifier that will be used as the column name for
            the forcing function trajectory in the returned data frame.

        Returns
        =======
        pandas.DataFrame
            A data frame indexed by time with all of the coordinates and
            measurements as columns.

        """
        t = self._calc_times(final_time, initial_time, sample_rate)

        typ = self._solution_type()

        x0, v0 = self._initial_conditions()

        m, c, k = self._canonical_coefficients()

        Fo = amplitude
        w = frequency
        fo = Fo / m

        if typ == 'no_damping':
            wn = self._natural_frequency(m, k)
            if math.isclose(w, wn):
                X = fo / 2 / w
                xss = X * t * np.sin(w*t)
                vss = X * w * t * np.cos(w*t) + X * np.sin(w*t)
                ass = (-X * w**2 * t * np.sin(w*t) + X * w * np.cos(w*t) + X *
                       w * np.cos(w*t))
                x = v0 / w * np.sin(w*t) + x0 * np.cos(w*t)
                v = v0 / w * w * np.cos(w*t) - x0 * w * np.sin(w*t)
                a = -v0 / w * w**2 * np.sin(w*t) - x0 * w**2 * np.cos(w*t)
            else:
                # steady state solution (particular solution)
                X = fo / (wn**2 - w**2)
                xss = X * np.cos(w * t)
                vss = -X * w * np.sin(w * t)
                ass = -X * w**2 * np.cos(w * t)
                # transient solution (homogenous solution)
                A1 = v0 / wn  # sin
                A2 = x0 - fo / (wn**2 - w**2)  # cos
                A = np.sqrt(A1**2 + A2**2)
                phi = np.arctan2(A2, A1)
                x = A * np.sin(wn * t + phi)
                v = A * wn * np.cos(wn * t + phi)
                a = -A * wn**2 * np.sin(wn * t + phi)
        elif typ == 'underdamped':
            wn, z, wd = self._normalized_form(m, c, k)

            theta = np.arctan2(2*z*wn*w, wn**2 - w**2)
            X = fo / np.sqrt((wn**2 - w**2)**2 + (2*z*wn*w)**2)

            xss = X * np.cos(w*t - theta)
            vss = -X * w * np.sin(w*t - theta)
            ass = -X * w**2 * np.cos(w*t - theta)

            phi = np.arctan2(wd * (x0 - X * np.cos(theta)),
                             v0 + (x0 - X * np.cos(theta)) *
                             z * wn - w * X * np.sin(theta))
            A = (x0 - X * np.cos(theta)) / np.sin(phi)
            x, v, a = self._damped_sinusoid(A, phi, t)
        else:
            raise ValueError('{} not yet supported.'.format(typ))

        self.result = self._state_traj_to_dataframe(t, x + xss, v + vss,
                                                    a + ass)
        if not col_name.isidentifier():
            msg = "'{}' is not a valid Python identifier."
            raise ValueError(msg.format(col_name))
        elif col_name in self.result.columns:
            raise ValueError('{} already taken.'.format(col_name))
        else:
            self.result[col_name] = amplitude * np.cos(frequency * t)

        return self.result

    def frequency_response(self, frequencies, amplitude):
        """Returns the amplitude and phase shift for simple sinusoidal forcing
        of the system. The first holds the plot of the coordinate's amplitude
        as a function of forcing frequency and the second holds a plot of the
        coordinate's phase shift with respect to the forcing function.

        Parameters
        ==========
        frequencies : array_like, shape(n,)
        amplitude : float
            The value of the forcing amplitude.

        Returns
        =======
        amp_curve : ndarray, shape(n,)
            The amplitude values of the coordinate at different frequencies.
        phase_curve : ndarray, shape(n,)
            The phase shift values in radians of the coordinate relative to the
            forcing.

        """
        m, c, k = self._canonical_coefficients()
        wn, z, wd = self._normalized_form(m, c, k)

        fo = amplitude / m
        w = np.asarray(frequencies)

        amp_curve = fo / np.sqrt((wn**2 - w**2)**2 + (2*z*wn*w)**2)
        phase_curve = np.arctan2(2*z*wn*w, wn**2 - w**2)

        return amp_curve, phase_curve

    def frequency_response_plot(self, amplitude, log=False, axes=None):
        """Returns an array of two matplotlib axes. The first holds the plot of
        the coordinate's amplitude as a function of forcing frequency and the
        second holds a plot of the coordinate's phase shift with respect to the
        forcing function.

        Parameters
        ==========
        amplitude : float
            The value of the forcing amplitude.
        log : boolean, optional
            If True, the amplitude will be plotted on a semi-log Y plot.

        """
        m, c, k = self._canonical_coefficients()
        wn, z, wd = self._normalized_form(m, c, k)

        w = np.linspace(0.0, 5 * wn, num=200)

        amp_curve, phase_curve = self.frequency_response(w, amplitude)

        if axes is None:
            fig, axes = plt.subplots(2, 1, sharex=True)
        axes[0].axvline(wn, color='black')
        if log:
            axes[0].semilogy(w, amp_curve)
        else:
            axes[0].plot(w, amp_curve)
        axes[0].set_ylabel('Coordinate Amplitude')
        axes[1].axvline(wn, color='black')
        axes[1].plot(w, np.rad2deg(phase_curve))
        axes[1].set_ylabel('Phase Shift [deg]')
        axes[1].set_xlabel('Forcing Frequency, $\omega$, [rad/s]')

        return axes


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
            the length of the edge of the book which is tangent to the cup's
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
            the angular rate of the book with respect to the gravity vector

    """

    def __init__(self):

        super(BookOnCupSystem, self).__init__()

        self.constants['thickness'] = 0.029  # m
        self.constants['length'] = 0.238  # m
        self.constants['radius'] = 0.042  # m
        self.constants['mass'] = 1.058  # kg
        self.coordinates['book_angle'] = 0.0  # rad
        self.speeds['book_angle_vel'] = 0.0  # rad/s

        def coeffs(thickness, length, radius):
            """Students will write this function themselves and pass it into
            the class via self.add_coef_func() when they get to modeling."""
            g = 9.81
            m = thickness**2 / 3 + length**2 / 12
            c = 0.0
            k = g * radius - g * thickness / 2
            return m, c, k

        self.canonical_coeffs_func = coeffs


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
            dissipation from things like air resistance, slip, etc.
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

        def coeffs(rotational_inertia, torsional_damping, torsional_stiffness):
            return rotational_inertia, torsional_damping, torsional_stiffness

        self.canonical_coeffs_func = coeffs


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

        def coeffs(pendulum_mass, inertia_about_joint, joint_to_mass_center,
                   acc_due_to_gravity):
            m = pendulum_mass
            i = inertia_about_joint
            l = joint_to_mass_center
            g = acc_due_to_gravity

            return i, 0.0, m * g * l

        self.canonical_coeffs_func = coeffs


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

        def coeffs(pendulum_mass, pendulum_length, acc_due_to_gravity):
            m = pendulum_mass
            l = pendulum_length
            g = acc_due_to_gravity

            return m * l**2, 0.0, m * g * l

        self.canonical_coeffs_func = coeffs


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

        def coeffs(bob_mass, bob_radius, rod_mass, rod_length, viscous_damping,
                   acc_due_to_gravity):

            Irod_O = rod_mass * rod_length**2 / 3
            Ibob_P = bob_mass * bob_radius**2 / 2
            Ibob_O = Ibob_P + bob_mass * rod_length**2

            I = Irod_O + Ibob_O
            C = viscous_damping * rod_length**2
            K = acc_due_to_gravity * rod_length * (bob_mass + rod_mass / 2.0)

            return I, C, K

        self.canonical_coeffs_func = coeffs


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

        self.coordinates['position'] = 0.0
        self.speeds['velocity'] = 0.0

        def coeffs(mass, damping, stiffness):
            return mass, damping, stiffness

        self.canonical_coeffs_func = coeffs


class BaseExcitationSystem(SingleDoFLinearSystem):
    """This system represents a mass connected to a moving massless base via a
    spring and damper in parallel. The motion of the mass is subject to viscous
    damping. The system is described by:

    Attributes
    ==========
    constants
        mass, m [kg]
            The suspended mass.
        damping, c [kg / s]
            The viscous linear damping coefficient which represents any energy
            dissipation from things like air resistance, friction, etc.
        stiffness, k [N / m]
            The linear elastic stiffness of the spring.
    coordinates
        position, x [m]
            The absolute position of the mass.
    speeds
        velocity, x_dot [m / s]
            The absolute velocity of the mass.

    """

    def __init__(self):

        super(BaseExcitationSystem, self).__init__()

        self.constants['mass'] = 1.0  # m
        self.constants['damping'] = 0.0  # kg/s
        self.constants['stiffness'] = 100  # N/m

        self.coordinates['position'] = 0.0
        self.speeds['velocity'] = 0.0

        def coeffs(mass, damping, stiffness):
            return mass, damping, stiffness

        self.canonical_coeffs_func = coeffs

    def sinusoidal_base_displacing_response(self, amplitude, frequency,
                                            final_time, initial_time=0.0,
                                            sample_rate=100,
                                            force_col_name='forcing_function',
                                            displace_col_name='displacing_function'):
        """Returns the trajectory of the system's coordinates, speeds,
        accelerations, and measurements if a sinusoidal displacement function
        described by:

        y(t) = Y * sin(ω*t)

        is specified for the movement of the base in the direction of the
        system's coordinate.

        Parameters
        ==========
        amplitude : float
            The amplitude of the displacement function, Y, in meters.
        frequency : float
            The frequency, ω, in radians per second of the sinusoidal
            displacement.
        final_time : float
            A value of time in seconds corresponding to the end of the
            simulation.
        initial_time : float, optional
            A value of time in seconds corresponding to the start of the
            simulation.
        sample_rate : integer, optional
            The sample rate of the simulation in Hertz (samples per second).
            The time values will be reported at the initial time and final
            time, i.e. inclusive, along with times space equally based on the
            sample rate.
        force_col_name : string, optional
            A valid Python identifier that will be used as the column name for
            the forcing function trajectory in the returned data frame.
        displace_col_name : string, optional
            A valid Python identifier that will be used as the column name for
            the forcing function trajectory in the returned data frame.

        Returns
        =======
        pandas.DataFrame
            A data frame indexed by time with all of the coordinates and
            measurements as columns.

        """
        m, c, k = self._canonical_coefficients()

        a0 = 0.0
        a1 = c * amplitude * frequency
        b1 = k * amplitude

        self.periodic_forcing_response(a0, a1, b1, frequency, final_time,
                                       initial_time=initial_time,
                                       sample_rate=sample_rate,
                                       col_name=force_col_name)

        try:
            col_name = self._displace_col_name
        except AttributeError:
            col_name = displace_col_name
        else:
            # if not the default warn
            if displace_col_name != 'displacing_function':
                msg = 'displace_col_name set to {}'
                warnings.warn(msg.format(displace_col_name))

        if col_name.isidentifier():
            self.result[col_name] = amplitude * np.sin(frequency *
                                                       self.result.index)
        else:
            msg = "'{}' is not a valid Python identifier."
            raise ValueError(msg.format(col_name))

        return self.result

    def periodic_base_displacing_response(self, twice_avg, cos_coeffs,
                                          sin_coeffs, frequency, final_time,
                                          initial_time=0.0, sample_rate=100,
                                          force_col_name='forcing_function',
                                          displace_col_name='displacing_function'):
        """Returns the trajectory of the system's coordinates, speeds,
        accelerations, and measurements if a periodic function defined by a
        Fourier series is applied as displacement of the base in the same
        direction as the system's coordinate. The displacing function is
        defined as::

                             N
            y(t)  = a0 / 2 + ∑ (an * cos(n*ω*t) + bn * sin(n*ω*t))
                            n=1

        Where a0, a1...an, and b1...bn are the Fourier coefficients. If N=∞
        then the Fourier series can describe any periodic function with a
        period (2*π)/ω.

        Parameters
        ==========
        twice_avg : float
            Twice the average value over one cycle, a0.
        cos_coeffs : float or sequence of floats
            The N cosine Fourier coefficients: a1, ..., aN.
        sin_coeffs : float or sequence of floats
            The N sine Fourier coefficients: b1, ..., bN.
        frequency : float
            The frequency, ω, in radians per second corresponding to one full
            cycle of the function.
        final_time : float
            A value of time in seconds corresponding to the end of the
            simulation.
        initial_time : float, optional
            A value of time in seconds corresponding to the start of the
            simulation.
        sample_rate : integer, optional
            The sample rate of the simulation in Hertz (samples per second).
            The time values will be reported at the initial time and final
            time, i.e. inclusive, along with times space equally based on the
            sample rate.
        force_col_name : string, optional
            A valid Python identifier that will be used as the column name for
            the forcing function trajectory in the returned data frame.
        displace_col_name : string, optional
            A valid Python identifier that will be used as the column name for
            the forcing function trajectory in the returned data frame.

        Returns
        =======
        pandas.DataFrame
            A data frame indexed by time with all of the coordinates, speeds,
            measurements, and forcing/displacing functions as columns.

        """
        t = self._calc_times(final_time, initial_time, sample_rate)

        m, c, k = self._canonical_coefficients()

        # shape(N,)
        cos_coeffs = np.atleast_1d(cos_coeffs)
        sin_coeffs = np.atleast_1d(sin_coeffs)

        N = len(cos_coeffs)
        # shape(N,)
        n = np.arange(1, N+1)

        a0 = k * twice_avg
        # shape(N,)
        an = k * cos_coeffs + c * sin_coeffs * n * frequency
        bn = k * sin_coeffs - c * cos_coeffs * n * frequency

        self.periodic_forcing_response(a0, an, bn, frequency, final_time,
                                       initial_time=initial_time,
                                       sample_rate=sample_rate,
                                       col_name=force_col_name)
        # shape(N, 1)
        n = np.arange(1, N+1)[:, np.newaxis]
        cos_coeffs = np.atleast_2d(cos_coeffs).T
        sin_coeffs = np.atleast_2d(sin_coeffs).T
        ycn = cos_coeffs * np.cos(n * frequency * t)
        ysn = sin_coeffs * np.sin(n * frequency * t)
        y = twice_avg / 2 + np.sum(ycn + ysn, axis=0)

        try:
            col_name = self._displace_col_name
        except AttributeError:
            col_name = displace_col_name
        else:
            # if not the default warn
            if displace_col_name != 'displacing_function':
                msg = 'displace_col_name set to {}'
                warnings.warn(msg.format(displace_col_name))

        if col_name.isidentifier():
            self.result[col_name] = y
        else:
            msg = "'{}' is not a valid Python identifier."
            raise ValueError(msg.format(col_name))

        return self.result


class SimpleQuarterCarSystem(BaseExcitationSystem):
    """This system represents a mass connected to a moving massless base via a
    spring and damper in parallel. The motion of the mass is subject to viscous
    damping. The system is described by:

    Attributes
    ==========
    constants
        mass, m [kg]
            The suspended mass.
        damping, c [kg / s]
            The viscous linear damping coefficient which represents any energy
            dissipation from things like air resistance, friction, etc.
        stiffness, k [N / m]
            The linear elastic stiffness of the spring.
    coordinates
        position, x [m]
            The absolute position of the mass.
    speeds
        velocity, x_dot [m / s]
            The absolute velocity of the mass.

    """
    _displace_col_name = 'road_height'

    def __init__(self):

        # NOTE : Don't call init on this super class, but the super super
        # class.
        super(BaseExcitationSystem, self).__init__()

        self.constants['sprung_mass'] = 1007  # kg
        self.constants['suspension_damping'] = 20E2  # kg/s
        self.constants['suspension_stiffness'] = 4E4  # N/m
        self.constants['travel_speed'] = 7.5  # m/s

        self.coordinates['car_vertical_position'] = -0.05  # m
        self.speeds['car_vertical_velocity'] = 0.0  # m/s

        xeq = 0.1  # m
        view_width = 4  # m
        # a rectangle will represent the car
        rect_width = 1.0  # m
        rect_height = rect_width / 6  # m

        def plot_config(car_vertical_position, time):

            fig, ax = plt.subplots(1, 1)

            ax.set_ylim((-0.1, 0.6))
            ax.set_ylabel('Height [m]')

            lat_pos = 0  # m

            lat = np.linspace(lat_pos - view_width / 2,
                              lat_pos + view_width / 2,
                              num=100)

            ax.set_xlim((lat[0], lat[-1]))

            rect = Rectangle(
                            (-rect_width / 2, xeq + car_vertical_position),  # (x,y)
                            rect_width,  # width
                            rect_height,  # height
                            )

            ax.add_patch(rect)

            time_txt = ax.text(lat_pos, 0.5, 'Time: {:1.1f} s'.format(time))

            # NOTE: Just plot a flat road for now because there may be no
            # results available
            road = ax.plot(lat, np.zeros_like(lat), color='black')[0]

            suspension = ax.plot([lat_pos, lat_pos],
                                 [0.0, xeq + car_vertical_position],
                                 linewidth='4', marker='o', color='yellow')[0]
            #force_vec = ax.plot([lat_pos, lat_pos],
                                #[xeq + car_vertical_position + rect_height / 2,
                                 #xeq + car_vertical_position + rect_height / 2 + 0.2],
                                #'r', linewidth=4)[0]

            return fig, ax, rect, road, suspension, time_txt

        self.config_plot_func = plot_config

        def plot_update(travel_speed, car_vertical_position,
                        time, time__hist, time__futr,
                        road_height, road_height__hist, road_height__futr,
                        ax, rect, road, suspension, time_txt):

            # v is a function of forcing freq
            lat_pos = travel_speed * time

            time_txt.set_text('Time: {:1.1f} s'.format(time))
            time_txt.set_position((lat_pos, 0.5))

            ax.set_xlim((lat_pos - view_width / 2, lat_pos + view_width / 2))

            rect.set_xy([lat_pos - rect_width / 2, xeq + car_vertical_position])

            lat = travel_speed * np.hstack((time__hist, time__futr))

            road.set_xdata(lat)
            road.set_ydata(np.hstack((road_height__hist, road_height__futr)))

            suspension.set_xdata([lat_pos, lat_pos])
            suspension.set_ydata([road_height, xeq + car_vertical_position])

            #force_vec.set_xdata([lat_pos, lat_pos])
            #force_vec.set_ydata([xeq + x[i] + rect_height / 2,
                                #xeq + x[i] + rect_height / 2 + f[i] / k])

        self.config_plot_update_func = plot_update

        def coeffs(sprung_mass, suspension_damping, suspension_stiffness):
            return sprung_mass, suspension_damping, suspension_stiffness

        self.canonical_coeffs_func = coeffs
