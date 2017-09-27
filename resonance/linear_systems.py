import math
import collections
from inspect import getargspec

import numpy as np
import pandas as pd
import matplotlib.animation as animation


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
    coordinates : _ParametersDict
        A custom dictionary that contains parameters that do vary with time.
    measurements : _MeasurementsDict
        A custom dictionary that contains parameters that are functions of the
        constants and coordinates.
    config_plot_func : function
        A function that generates a matplotlib plot that uses the instantaneous
        values of the constants, coordinates, and measurements.
    config_plot_update_func : function
        A function that updates the configuration plot that uses the time
        series values of the constants, coordinates, and measurements. Defines
        an matplotlib animation frame.

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

    """

    _time_var_name = 'time'

    # TODO : Only allow a single coordinate to be set on a 1 DoF system.

    def __init__(self):

        # TODO : Allow pars, coords, and meas to be set on intialization.

        self._constants = _ParametersDict({})
        self._coordinates = _ParametersDict({})
        self._measurements = _MeasurementsDict({})

        self._measurements._constants = self._constants
        self._measurements._coordinates = self._coordinates
        self._measurements._funcs = {}

        self._config_plot_func = None
        self._config_plot_update_func = None

    @property
    def constants(self):
        """A dictionary containing the all of the system's constants."""
        return self._constants

    @constants.setter
    def constants(self, val):
        msg = ('It is not allowed to replace the entire constants dictionary, '
               'add or delete constants one by one.')
        raise ValueError(msg)

    @property
    def coordinates(self):
        """A dictionary containing the all of the system's coordinates."""
        return self._coordinates

    @coordinates.setter
    def coordinates(self, val):
        msg = ('It is not allowed to replace the entire coordinates '
               'dictionary, add or delete coordinates one by one.')
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
        return self._config_plot_func

    @config_plot_func.setter
    def config_plot_func(self, func):
        """The configuration plot function arguments should be any of the
        system's constants, coordinates, measurements, or 'time'. No other
        arguments are valid. The function has to return the matplotlib figure
        as the first item but can be followed by any number of mutable
        matplotlib objects that you may want to change during an animation."""
        self._config_plot_func = func

    @property
    def config_plot_update_func(self):
        return self._config_plot_update_func

    @config_plot_update_func.setter
    def config_plot_update_func(self, func):
        """The configuration plot update function arguments should be any of
        the system's constants, coordinates, measurements, or 'time'. No other
        arguments are valid. Nothing need be returned from the function."""
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
            with any names in the parameters or coordinates dictionary.
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
                      list(self.coordinates.keys())):
            msg = ('{} is already used as a constant or coordinate name. '
                   'Choose something different.')
            raise ValueError(msg.format(name))

        self.measurements._funcs[name] = func
        dict.__setitem__(self.measurements, name,
                         self.measurements._compute_value(name))

    @staticmethod
    def _natural_frequency(mass, stiffness):
        return np.sqrt(stiffness / mass)

    @staticmethod
    def _damping_ratio(mass, damping, natural_frequency):
        return damping / 2.0 / mass / natural_frequency

    @staticmethod
    def _damped_natural_frequency(natural_frequency, damping_ratio):
        return natural_frequency * np.sqrt(1.0 - damping_ratio**2)

    def _solution_func(self):

        m, c, k = self._canonical_coefficients()

        if math.isclose(c, 0.0):
            sol_func = self._no_damping_solution
        else:  # damping, so check zeta
            omega_n = self._natural_frequency(m, k)
            zeta = self._damping_ratio(m, c, omega_n)
            if zeta < 1.0:
                sol_func = self._underdamped_solution
            elif zeta > 1.0:
                sol_func = self._overdamped_solution
            elif math.isclose(zeta, 0.0):
                sol_func = self._critically_damped_solution
            else:
                sol_func = None

        return sol_func

    def _no_damping_solution(self, time):

        t = time

        m, c, k = self._canonical_coefficients()

        wn = self._natural_frequency(m, k)

        x0 = list(self.coordinates.values())[0]
        v0 = 0.0

        c1 = v0 / wn
        c2 = x0

        return c1 * np.sin(wn * t) + c2 * np.cos(wn * t)

    def _underdamped_solution(self, time):

        t = time

        m, c, k = self._canonical_coefficients()

        wn = self._natural_frequency(m, k)
        z = self._damping_ratio(m, c, wn)
        wd = self._damped_natural_frequency(wn, z)

        x0 = list(self.coordinates.values())[0]
        v0 = 0.0

        A = np.sqrt(((v0 + z * wn * x0)**2 + (x0 * wd)**2) / wd**2)
        phi = np.atan(x0 * wd / (v0 + z * wn * x0))

        return A * np.exp(-z * wn * t) * np.sin(wd * t + phi)

    def _overdamped_solution(self, time):

        t = time

        m, c, k = self._canonical_coefficients()

        wn = self._natural_frequency(m, k)
        z = self._damping_ratio(m, c, wn)

        x0 = list(self.coordinates.values())[0]
        v0 = 0.0

        a1 = ((-v0 + (-z + np.sqrt(z**2 - 1)) * wn * x0) / 2 / wn /
              np.sqrt(z**2 - 1))
        a2 = ((v0 + (z + np.sqrt(z**2 - 1)) * wn * x0) / 2 / wn /
              np.sqrt(z**2 - 1))

        decay = wn * np.sqrt(z**2 - 1) * t
        return np.exp(-z * wn * t) * (a1 * np.exp(-decay) + a2 * np.exp(decay))

    def _critically_damped_solution(self, time):

        t = time

        m, c, k = self._canonical_coefficients()

        wn = self._natural_frequency(m, k)

        x0 = list(self.coordinates.values())[0]
        v0 = 0.0

        a1 = x0
        a2 = v0 + wn * x0

        return (a1 + a2 * t) * np.exp(-wn * t)

    def period(self):
        """Returns the period of oscillation of the coordinate."""
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
            A data frame index by time with all of the coordinates and
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

        coordinate_traj = sol_func(times)

        coord_name = list(self.coordinates.keys())[0]

        df = pd.DataFrame({coord_name: coordinate_traj}, index=times)
        df.index.name = self._time_var_name

        # TODO : Need a way to compute the measurements as array values based
        # on the coordinate changing at each time but the current method of
        # letting the measurement be computed by the stored coordinate scalar
        # is a bit problematic.
        for k, v in self.measurements.items():
            vals = np.zeros_like(times)
            x0 = list(self.coordinates.values())[0]
            for i, xi in enumerate(coordinate_traj):
                self.coordinates[coord_name] = xi
                vals[i] = self.measurements[k]
            self.coordinates[coord_name] = x0
            df[k] = vals

        self.result = df

        return df

    def plot_configuration(self):
        """Returns a matplotlib figure generated by the function assigned to
        the config_plot_func attribute. You may need to call
        ``matplotlib.pyplot.show()`` to display the figure.

        Returns
        =======
        fig : matplotlib.figure.Figure

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

        if self.config_plot_func is None:
            msg = 'No ploting function has been assigned to config_plot_func.'
            raise ValueError(msg)

        if self.config_plot_update_func is None:
            msg = ('No ploting update function has been assigned to '
                   'config_plot_update_func.')
            raise ValueError(msg)

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

    """

    def __init__(self):

        super(BookOnCupSystem, self).__init__()

        self.constants['thickness'] = 0.029  # m
        self.constants['length'] = 0.238  # m
        self.constants['radius'] = 0.042  # m
        self.constants['mass'] = 1.058  # kg
        self.coordinates['book_angle'] = 0.0  # rad

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
