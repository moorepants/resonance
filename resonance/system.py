import collections as _collections
from inspect import getargspec

import numpy as np
import matplotlib.animation as animation
import pandas as pd


class _ParametersDict(_collections.MutableMapping, dict):
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


class _StatesDict(_collections.OrderedDict):
    def __setitem__(self, key, value, allow=False):
        msg = ("You can't set the values on the states dictionary, set them on "
               "the coordinates or speeds dictionaries instead.")
        if allow:
            _collections.OrderedDict.__setitem__(self, key, value)
        else:
            raise ValueError(msg)


class _CoordinatesDict(_collections.MutableMapping, dict):
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


class _MeasurementsDict(_collections.MutableMapping, dict):

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


class System(object):
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
    def states(self):
        """An ordered dictionary containing the system's state variables and
        values. The coordinates are always ordered before the speeds and the
        individual order of the values depends on the order they were added to
        coordinates and speeds."""
        states = _StatesDict({})
        for k, v in self.coordinates.items():
            states.__setitem__(k, v, allow=True)
        for k, v in self.speeds.items():
            states.__setitem__(k, v, allow=True)
        return states

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
        msg = '{} is not in constants, coordinates, speeds, or measurements.'
        try:
            v = self.constants[par_name]
        except KeyError:
            try:
                v = self.coordinates[par_name]
            except KeyError:
                try:
                    v = self.speeds[par_name]
                except KeyError:
                    try:
                        v = self.measurements[par_name]
                    except KeyError:
                        raise KeyError(msg.format(par_name))
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

    def _state_traj_to_dataframe(self, times, pos, vel, acc):

        coord_name = list(self.coordinates.keys())[0]
        speed_name = coord_name + self._vel_append

        # TODO : What if they added a coordinate with the vel or acc names?
        df = pd.DataFrame({coord_name: pos,
                           speed_name: vel,
                           coord_name + self._acc_append: acc},
                          index=times)
        df.index.name = self._time_var_name

        # TODO : Allow acc to be used in measurements.
        # TODO : Need a way to compute the measurements as array values based
        # on the coordinate changing at each time but the current method of
        # letting the measurement be computed by the stored coordinate scalar
        # is a bit problematic.
        for k, v in self.measurements.items():
            vals = np.zeros_like(times)
            x0 = list(self.coordinates.values())[0]
            v0 = list(self.speeds.values())[0]
            for i, (xi, vi) in enumerate(zip(pos, vel)):
                self.coordinates[coord_name] = xi
                self.speeds[speed_name] = vi
                vals[i] = self.measurements[k]
            self.coordinates[coord_name] = x0
            self.speeds[speed_name] = v0
            df[k] = vals

        return df

    def free_response(self, final_time, initial_time=0.0, sample_rate=100,
                      **kwargs):
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

        pos, vel, acc = self._generate_state_trajectories(times, **kwargs)

        self.result = self._state_traj_to_dataframe(times, pos, vel, acc)

        return self.result

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
                elif k == 'time__hist':
                    args.append(0.0)  # static config defaults to t=0.0
                elif k.endswith('__hist'):
                    args.append(self._get_par_vals(k[:-6]))
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
                elif k == 'time__hist':
                    args.append(self.result[:time].index)
                elif k.endswith('__hist'):
                    args.append(self.result[k[:-6]][:time])
                else:
                    try:
                        args.append(row[k])
                    except KeyError:
                        # requires these to be in the same order
                        args.append(pop_list.pop(0))
            self.config_plot_update_func(*args)

        return animation.FuncAnimation(fig, gen_frame, fargs=(objs_to_modify, ),
                                       frames=self.result.iterrows(), **kwargs)
