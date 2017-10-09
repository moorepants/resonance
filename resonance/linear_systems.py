import math
import collections
from inspect import getargspec
import warnings

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib.patches import Circle


class _ParametersDict(collections.MutableMapping, dict):
    """A custom dictionary for storing constants and coordinates."""

    def __getitem__(self, key):
        return dict.__getitem__(self, key)

    def __setitem__(self, key, value):
        # TODO : It would be nice to check to ensure that these names do not
        # clash with names in measurements or coordinates. Not quite sure how
        # one would do that.
        if not key.isidentifier():
            msg = ('{} is not a valid parameter name. '
                   'It must be a valid Python variable name.')
            raise ValueError(msg.format(key))
        elif key.lower() == 'time':
            msg = ('{} is a reserved parameter name. '
                   'Choose something different.')
            raise ValueError(msg.format(key))
        else:
            dict.__setitem__(self, key, value)

    def __delitem__(self, key):
        dict.__delitem__(self, key)

    def __iter__(self):
        return dict.__iter__(self)

    def __len__(self):
        return dict.__len__(self)

    def __contains__(self, x):
        return dict.__contains__(self, x)


class _CoordinatesDict(collections.MutableMapping, dict):
    """A custom dictionary for storing coordinates and speeds."""

    def __getitem__(self, key):
        return dict.__getitem__(self, key)

    def __setitem__(self, key, value):
        if len(self.keys()) == 1 and key != list(self.keys())[0]:
            msg = ("There is already a coordinate set for this system, only "
                   "one coordinate is permitted. Use del to remove the "
                   "coordinate if you'd like to add a new one.")
            raise ValueError(msg)
        elif not key.isidentifier():
            msg = ('{} is not a valid coordinate or speed name. '
                   'It must be a valid Python variable name.')
            raise ValueError(msg.format(key))
        elif key.lower() == 'time':
            msg = ('{} is a reserved parameter name. '
                   'Choose something different.')
            raise ValueError(msg.format(key))
        else:
            dict.__setitem__(self, key, value)

    def __delitem__(self, key):
        dict.__delitem__(self, key)

    def __iter__(self):
        return dict.__iter__(self)

    def __len__(self):
        return dict.__len__(self)

    def __contains__(self, x):
        return dict.__contains__(self, x)


class _MeasurementsDict(collections.MutableMapping, dict):

    def _compute_value(self, key):

        func = self._funcs[key]

        def get_par(k):
            try:
                v = self._constants[k]
            except KeyError:
                try:
                    v = self._coordinates[k]
                except KeyError:
                    try:
                        v = self._speeds[k]
                    except KeyError:
                        v = self[k]
            return v

        # TODO : getargspec is deprecated, supposedly signature cna do the same
        # thing but the args are in a dictionary and it isn't clear to me they
        # are ordered.
        args = [get_par(k) for k in getargspec(func).args]
        return func(*args)

    def __getitem__(self, key):
        return self._compute_value(key)

    def __setitem__(self, key, value):
        msg = ('It is not possible to assign a value to a measurement, '
               'add a measurement function instead.')
        raise ValueError(msg)

    def __delitem__(self, key):
        # TODO : Should also remove from self._funcs
        dict.__delitem__(self, key)

    def __iter__(self):
        return dict.__iter__(self)

    def __len__(self):
        return dict.__len__(self)

    def __contains__(self, x):
        return dict.__contains__(self, x)


class SingleDoFLinearSystem(object):
    """This is the base system for any linear single degree of freedom system.
    It can be subclassed to make a custom system or instantiated and the
    attributes and methods added dynamically.

    Attributes
    ==========
    constants : _ParametersDict
        A custom dictionary that contains parameters that do not vary with
        time.
    coordinates : _CoordinatesDict
        A custom dictionary that contains the generalized coordinate which
        varies with time.
    speeds : _CoordinatesDict
        A custom dictionary that contains the generalized speed which varies
        with time.
    measurements : _MeasurementsDict
        A custom dictionary that contains parameters that are functions of the
        constants, coordinates, and other measurements.
    config_plot_func : function
        A function that generates a matplotlib plot that uses the instantaneous
        values of the constants, coordinates, and measurements.
    config_plot_update_func : function
        A function that updates the configuration plot that uses the time
        series values of the constants, coordinates, and measurements. Defines
        a matplotlib animation frame.

    Methods
    =======
    add_measurement
        Used to dynamically add functions to compute measurement values.
    free_response
        Simulates the system and returns a time series of the coordinates and
        measurments.
    plot_configuration
        Generates the plot defined by ``config_plot_func``.
    animate_configutation
        Generates the animation defined by ``config_plot_func`` and
        ``config_plot_update_func``.
    period
        Returns the damped natural period of the system.

    """

    _time_var_name = 'time'
    _vel_append = '_vel'
    _acc_append = '_acc'

    def __init__(self):

        # TODO : Allow constants, coords, and meas to be set on intialization.

        self._constants = _ParametersDict({})
        self._coordinates = _CoordinatesDict({})
        self._speeds = _CoordinatesDict({})
        self._measurements = _MeasurementsDict({})

        self._measurements._constants = self._constants
        self._measurements._coordinates = self._coordinates
        self._measurements._speeds = self._speeds
        self._measurements._funcs = {}

        self._config_plot_func = None
        self._config_plot_update_func = None

    @property
    def constants(self):
        """A dictionary containing the all of the system's constants.

        Examples
        ========
        >>> sys = SingleDoFLinearSystem()
        >>> sys.constants
        {}
        >>> sys.constants['mass'] = 5.0
        >>> sys.constants
        {'mass': 5.0}
        >>> del sys.constants['mass']
        >>> sys.constants
        {}
        >>> sys.constants['length'] = 10.0
        >>> sys.constants
        {'length': 10.0}


        """
        return self._constants

    @constants.setter
    def constants(self, val):
        msg = ('It is not allowed to replace the entire constants dictionary, '
               'add or delete constants one by one.')
        raise ValueError(msg)

    @property
    def coordinates(self):
        """A dictionary containing the system's generalized coordinate."""
        return self._coordinates

    @coordinates.setter
    def coordinates(self, val):
        msg = ('It is not allowed to replace the entire coordinates '
               'dictionary, add or delete coordinates one by one.')
        raise ValueError(msg)

    @property
    def speeds(self):
        """A dictionary containing the system's generalized speed."""
        return self._speeds

    @speeds.setter
    def speeds(self, val):
        msg = ('It is not allowed to replace the entire speeds '
               'dictionary, add or delete speeds one by one.')
        raise ValueError(msg)

    @property
    def measurements(self):
        """A dictionary containing the all of the system's measurements."""
        return self._measurements

    @measurements.setter
    def measurements(self, val):
        msg = ('It is not allowed to replace the entire measurements '
               'dictionary; add measurement functions using the '
               'add_measurement() method.')
        raise ValueError(msg)

    @property
    def config_plot_func(self):
        """The configuration plot function arguments should be any of the
        system's constants, coordinates, measurements, or 'time'. No other
        arguments are valid. The function has to return the matplotlib figure
        as the first item but can be followed by any number of mutable
        matplotlib objects that you may want to change during an animation.
        Refer to the matplotlib documentation for tips on creating figures.

        Examples
        ========
        >>> sys = SingleDoFLinearSystem()
        >>> sys.constants['radius'] = 5.0
        >>> sys.constants['center_y'] = 10.0
        >>> sys.coordinates['center_x'] = 0.0
        >>> def plot(radius, center_x, center_y, time):
        ...     fig, ax = plt.subplots(1, 1)
        ...     circle = Circle((center_x, center_y), radius=radius)
        ...     ax.add_patch(circle)
        ...     ax.set_title(time)
        ...     return fig, circle, ax
        ...
        >>> sys.config_plot_function = plot
        >>> sys.plot_configuration()

        """
        return self._config_plot_func

    @config_plot_func.setter
    def config_plot_func(self, func):
        self._config_plot_func = func

    @property
    def config_plot_update_func(self):
        """The configuration plot update function arguments should be any of
        the system's constants, coordinates, measurements, or 'time' in any
        order with the returned values from the ``config_plot_func`` as the
        last arguments in the exact order as in the configuration plot return
        statement. No other arguments are valid. Nothing need be returned from
        the function. See the matplotlib animation documentation for tips on
        creating these update functions.

        Examples
        ========
        >>> sys = SingleDoFLinearSystem()
        >>> sys.constants['radius'] = 5.0
        >>> sys.constants['center_y'] = 10.0
        >>> sys.coordinates['center_x'] = 0.0
        >>> def plot(radius, center_x, center_y, time):
        ...     fig, ax = plt.subplots(1, 1)
        ...     circle = Circle((center_x, center_y), radius=radius)
        ...     ax.add_patch(circle)
        ...     ax.set_title(time)
        ...     return fig, circle, ax
        ...
        >>> sys.config_plot_function = plot
        >>> def update(center_y, center_x, time, circle, ax):
        ...     # NOTE : that circle and ax have to be the last arguments and be
        ...     # in the same order as returned from plot()
        ...     circle.set_xy((center_x, center_y))
        ...     ax.set_title(time)
        ...     fig.canvas.draw()
        ...
        >>> sys.config_plot_update_func = update
        >>> sys.animate_configuration()


        """
        return self._config_plot_update_func

    @config_plot_update_func.setter
    def config_plot_update_func(self, func):
        self._config_plot_update_func = func

    def _get_par_vals(self, par_name):
        """Returns the value of any variable stored in the parameters,
        coordinates, or measurements dictionaries."""
        # TODO : Maybe just merging the three dicts is better? What to do about
        # name clashes in the dictionaries? Garbage in garbage out?
        # TODO : Raise useful error message if par_name not in any of the
        # dicts.
        try:
            v = self.constants[par_name]
        except KeyError:
            try:
                v = self.coordinates[par_name]
            except KeyError:
                v = self.measurements[par_name]
        return v

    def add_measurement(self, name, func):
        """Creates a new measurement entry in the measurements attribute that
        uses the provided function to compute the measurement.

        Parameters
        ==========
        name : string
            This must be a valid Python variable name and it should not clash
            with any names in the constants or coordinates dictionary.
        func : function
            This function must only have existing parameter, coordinate, or
            measurement names in the function signature. These can be a subset
            of the available choices and any order is permitted. The function
            must be able to operate on arrys, i.e. use NumPy vectorized
            functions inside. It should return a single variable, scalar or
            array, that gives the values of the measurement.


        Examples
        ========

        >>> import numpy as np
        >>> def f(par2, meas4, par1, coord5):
                  return par2 + meas4 + par1 + np.abs(coord5)
        >>> f(1.0, 2.0, 3.0, -4.0):
        10.0
        >>> f(1.0, 2.0, 3.0, np.array([1.0, 2.0, -4.0]))
        array([  7.,   8.,  10.])
        >>> sys.add_measurement('meas5', f)
        >>> sys.measurements['meas5']
        10.0

        """
        if name.lower() == 'time':
            msg = ('{} is a reserved parameter name. '
                   'Choose something different.')
            raise ValueError(msg.format(name))
        elif name in (list(self.constants.keys()) +
                      list(self.coordinates.keys()) +
                      list(self.speeds.keys())):
            msg = ('{} is already used as a constant or coordinate name. '
                   'Choose something different.')
            raise ValueError(msg.format(name))

        self.measurements._funcs[name] = func
        dict.__setitem__(self.measurements, name,
                         self.measurements._compute_value(name))

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

    def free_response(self, final_time, initial_time=0.0, sample_rate=100):
        """Returns a Pandas data frame with monotonic time values as the index
        and columns for each coordinate and measurement at the time value for
        that row. Note that this data frame is stored on the system as the
        vairable ``result`` until this method is called again, which will
        overwrite it.

        Parameters
        ==========
        final_time : float
            A value of time in seconds corresponding to the end of the
            simulation.
        initial_time : float, optional
            A value of time in seconds corresponding to the start of the
            simulation.
        sample_rate : integer, optional
            The sample rate of the simulation in Hertz (samples per second).
            The time values will be reported at the initial time and final
            time, along with times space equally based on the sample rate.

        Returns
        =======
        df : pandas.DataFrame
            A data frame indexed by time with all of the coordinates and
            measurements as columns.

        """
        # TODO : Should have the option to pass in unequally spaced monotonic
        # time arrays.

        if final_time < initial_time:
            raise ValueError('Final time must be greater than initial time.')

        delta = final_time - initial_time
        num = int(round(sample_rate * delta)) + 1
        times = np.linspace(initial_time, final_time, num=num)

        sol_func = self._solution_func()

        pos_traj, vel_traj, acc_traj = sol_func(times)

        coord_name = list(self.coordinates.keys())[0]
        speed_name = coord_name + self._vel_append

        # TODO : What if they added a coordinate with the vel or acc names?
        df = pd.DataFrame({coord_name: pos_traj,
                           speed_name: vel_traj,
                           coord_name + self._acc_append: acc_traj},
                          index=times)
        df.index.name = self._time_var_name

        # TODO : Allow vel and acc to be used in measurements.
        # TODO : Need a way to compute the measurements as array values based
        # on the coordinate changing at each time but the current method of
        # letting the measurement be computed by the stored coordinate scalar
        # is a bit problematic.
        for k, v in self.measurements.items():
            vals = np.zeros_like(times)
            x0 = list(self.coordinates.values())[0]
            v0 = list(self.speeds.values())[0]
            for i, (xi, vi) in enumerate(zip(pos_traj, vel_traj)):
                self.coordinates[coord_name] = xi
                self.speeds[speed_name] = vi
                vals[i] = self.measurements[k]
            self.coordinates[coord_name] = x0
            self.speeds[speed_name] = v0
            df[k] = vals
        self.result = df

        return df

    def plot_configuration(self):
        """Returns a matplotlib figure generated by the function assigned to
        the ``config_plot_func`` attribute. You may need to call
        ``matplotlib.pyplot.show()`` to display the figure.

        Returns
        =======
        fig : matplotlib.figure.Figure
            The first returned object is always a figure.
        *args : matplotlib objects
            Any matplotlib objects can be returned after the figure.

        """
        # TODO : Most plots in pandas, etc return the axes not the figure. I
        # think the parent figure can always be gotten from an axis.
        if self.config_plot_func is None:
            msg = 'No ploting function has been assigned to config_plot_func.'
            raise AttributeError(msg)
        else:
            args = []
            for k in getargspec(self.config_plot_func).args:
                if k == 'time':
                    args.append(0.0)  # static config defaults to t=0.0
                else:
                    args.append(self._get_par_vals(k))
            return self.config_plot_func(*args)

    def animate_configuration(self, **kwargs):
        """Returns a matplotlib animation function based on the configuration
        plot and the configuration plot update function."""

        if self.config_plot_update_func is None:
            msg = ('No ploting update function has been assigned to '
                   'config_plot_update_func.')
            raise AttributeError(msg)

        # TODO : Could be:
        # axes, *objs_to_modify = ..
        # try:
        #   fig = axes.figure
        # except AttributeError:
        #   fig = axes[0].figure
        try:
            fig, *objs_to_modify = self.plot_configuration()
        except TypeError:
            msg = ('The configuration plot function does not return any objects'
                   ' to modify in the update function.')
            raise ValueError(msg)

        def gen_frame(row_tuple, pop_list):
            time = row_tuple[0]
            row = row_tuple[1]
            # Don't mutate the orginal list.
            pop_list = pop_list.copy()
            args = []
            for k in getargspec(self.config_plot_update_func).args:
                if k == 'time':
                    args.append(time)
                else:
                    try:
                        args.append(row[k])
                    except KeyError:
                        # requires these to be in the same order
                        args.append(pop_list.pop(0))
            self.config_plot_update_func(*args)

        return animation.FuncAnimation(fig, gen_frame, fargs=(objs_to_modify, ),
                                       frames=self.result.iterrows(), **kwargs)


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
