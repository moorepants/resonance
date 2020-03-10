import sys
import collections as _collections
from collections import abc as _abc
from inspect import getfullargspec
import random

import numpy as np
import matplotlib.animation as animation
import pandas as pd

_FORBIDDEN_SUFFIXES = ['_acc', '__hist', '__futr']


class _ConstantsDict(_abc.MutableMapping, dict):
    """A custom dictionary for storing constants."""

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
        elif any([key.endswith(s) for s in _FORBIDDEN_SUFFIXES]):
            msg = ('{} are reserved suffixes. '
                   'Choose something different.')
            raise ValueError(msg.format(_FORBIDDEN_SUFFIXES))
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

    def __str__(self):
        pairs = []
        for k, v in self.items():
            pairs.append("'{}': {}".format(k, v))
        return '{' + ', '.join(pairs) + '}'

    def __repr__(self):
        return self.__str__()


class _CoordinatesDict(_StatesDict):
    """A custom dictionary for storing coordinates and speeds."""

    def __setitem__(self, key, value):
        # TODO : Shouldn't allow coordinates to be named with suffix _vel.
        if not key.isidentifier():
            msg = ('{} is not a valid coordinate or speed name. '
                   'It must be a valid Python variable name.')
            raise ValueError(msg.format(key))
        elif key.lower() == 'time':
            msg = ('{} is a reserved parameter name. '
                   'Choose something different.')
            raise ValueError(msg.format(key))
        elif any([key.endswith(s) for s in _FORBIDDEN_SUFFIXES]):
            msg = ('{} are reserved suffixes. '
                   'Choose something different.')
            raise ValueError(msg.format(_FORBIDDEN_SUFFIXES))
        else:
            _collections.OrderedDict.__setitem__(self, key, value)


class _SingleDoFCoordinatesDict(_CoordinatesDict):
    def __setitem__(self, key, value):
        if len(self.keys()) == 1 and key != list(self.keys())[0]:
            msg = ("There is already a coordinate set for this system, only "
                   "one coordinate is permitted. Use del to remove the "
                   "coordinate if you'd like to add a new one.")
            raise ValueError(msg)
        else:
            super(_SingleDoFCoordinatesDict, self).__setitem__(key, value)


class _MeasurementsDict(_abc.MutableMapping, dict):

    def _check_for_duplicate_keys(self):
        c_keys = list(self._constants.keys())
        g_keys = list(self._coordinates.keys())
        s_keys = list(self._speeds.keys())
        m_keys = list(self.keys())
        all_keys = c_keys + g_keys + s_keys + m_keys
        if len(set(all_keys)) < len(all_keys):
            dups = []
            sorted_keys = sorted(all_keys)
            for i, k in enumerate(sorted_keys[1:]):
                if k == sorted_keys[i - 1]:
                    dups.append(k)
            msg = ("{} are duplicate keys in your system's parameters. "
                   "Duplicates are not allowed.")
            raise KeyError(msg.format(dups))

    def _compute_value(self, key):

        func = self._funcs[key]

        def get_par(k):
            if k == 'time':
                v = self._time['t']
            elif k in self._constants:
                v = self._constants[k]
            elif k in self._coordinates:
                v = self._coordinates[k]
            elif k in self._speeds:
                v = self._speeds[k]
            elif k in self:
                v = self[k]
            else:
                msg = ("{} not in constants, coordinates, speeds, or "
                       "measurements.")
                raise KeyError(msg.format(k))
            return v

        # TODO : getfullargspec is deprecated, supposedly signature can do the
        # same thing but the args are in a dictionary and it isn't clear to me
        # they are ordered.
        args = [get_par(k) for k in getfullargspec(func).args]
        return func(*args)

    def __getitem__(self, key):
        # TODO : This call to _check_for_duplicate_keys() is a big CPU time
        # hog.
        self._check_for_duplicate_keys()
        return self._compute_value(key)

    def __setitem__(self, key, value):
        msg = ('It is not possible to assign a value to a measurement, '
               'add a measurement function instead.')
        raise ValueError(msg)

    def __delitem__(self, key):
        dict.__delitem__(self._funcs, key)
        dict.__delitem__(self, key)

    def __iter__(self):
        return dict.__iter__(self)

    def __len__(self):
        return dict.__len__(self)

    def __contains__(self, x):
        return dict.__contains__(self, x)


class System(object):
    """This is the abstract base class for any single or multi degree of
    freedom system. It can be sub-classed to make a custom system or the
    necessary methods can be added dynamically."""

    _time_var_name = 'time'
    _vel_append = '_vel'
    _acc_append = '_acc'

    def __init__(self):

        # TODO : Allow constants, coords, and meas to be set on initialization.

        self._time = {'t': 0.0}  # this needs to be a mutable object

        self._constants = _ConstantsDict({})
        self._coordinates = _CoordinatesDict({})
        self._speeds = _CoordinatesDict({})
        self._measurements = _MeasurementsDict({})

        self._measurements._constants = self._constants
        self._measurements._coordinates = self._coordinates
        self._measurements._speeds = self._speeds
        self._measurements._time = self._time
        self._measurements._funcs = {}

        self._config_plot_func = None
        self._config_plot_update_func = None

    def __str__(self):
        text = """\
System name: {name}

Configuration plot function defined: {plot}
Configuration update function defined: {update}

Constants
=========
{constants}

Coordinates
===========
{coordinates}

Speeds
======
{speeds}

Measurements
============
{measurements}"""

        ls = lambda d: '\n'.join('{} = {:1.5f}'.format(k, v) for k, v in d.items())

        speed_tmpl = '{} = d({})/dt = {:1.5f}'
        speed_strs = []
        for gs, gc in zip(self.speeds.keys(), self.coordinates.keys()):
            speed_strs.append(speed_tmpl.format(gs, gc, self.speeds[gs]))

        return text.format(
            name=self.__class__.__name__,
            plot='True' if self.config_plot_func is not None else 'False',
            update='True' if self.config_plot_update_func is not None else 'False',
            constants=ls(self.constants),
            coordinates=ls(self.coordinates),
            speeds='\n'.join(speed_strs),
            measurements=ls(self.measurements))

    def __repr__(self):
        return self.__str__()

    @property
    def constants(self):
        """A dictionary containing the all of the system's constants, i.e.
        parameters that do not vary with time.

        Examples
        ========
        >>> sys = System()
        >>> sys.constants
        {}
        >>> sys.constants['mass'] = 5.0
        >>> sys.constants
        {'mass': 5.0}
        >>> del sys.constants['mass']
        >>> sys.constants
        {}
        >>> sys.constants['mass'] = 5.0
        >>> sys.constants['length'] = 10.0
        >>> sys.constants
        {'mass': 5.0, 'length': 10.0}

        """
        return self._constants

    @constants.setter
    def constants(self, val):
        msg = ('It is not allowed to replace the entire constants dictionary, '
               'add or delete constants one by one.')
        raise ValueError(msg)

    @property
    def coordinates(self):
        """A dictionary containing the system's generalized coordinates, i.e.
        coordinate parameters that vary with time. These values will be used as
        initial conditions in simulations.

        Examples
        ========
        >>> sys = System()
        >>> sys.coordinates['angle'] = 0.0
        >>> sys.coordinates
        {'angle': 0.0}

        """
        return self._coordinates

    @coordinates.setter
    def coordinates(self, val):
        msg = ('It is not allowed to replace the entire coordinates '
               'dictionary, add or delete coordinates one by one.')
        raise ValueError(msg)

    @property
    def speeds(self):
        """A dictionary containing the system's generalized speeds, i.e. speed
        parameters that vary with time. These values will be used as initial
        conditions in simulations.

        Examples
        ========
        >>> sys = System()
        >>> sys.speeds['angle_vel'] = 0.0
        >>> sys.speeds
        {'angle_vel': 0.0}

        """
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
        coordinates and speeds.

        Examples
        ========
        >>> sys = System()
        >>> sys.coordinates['angle'] = 0.2
        >>> sys.speeds['angle_vel'] = 0.1
        >>> sys.states
        {'angle': 0.2, 'angle_vel': 0.1}
        >>> list(sys.states.keys())
        ['angle', 'angle_vel']
        >>> list(sys.states.values())
        [0.2, 0.1]

        """
        states = _StatesDict({})
        for k, v in self.coordinates.items():
            states.__setitem__(k, v, allow=True)
        for k, v in self.speeds.items():
            states.__setitem__(k, v, allow=True)
        return states

    @property
    def measurements(self):
        """A dictionary containing the all of the system's measurements, i.e.
        parameters that are functions of the constants, coordinates, speeds,
        and other measurements."""
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
        system's constants, coordinates, measurements, or ``'time'``. No other
        arguments are valid. The function has to return the matplotlib figure
        as the first item but can be followed by any number of mutable
        matplotlib objects that you may want to change during an animation.
        Refer to the matplotlib documentation for tips on creating figures.

        Examples
        ========

        >>> import matplotlib.pyplot as plt
        >>> from resonance.linear_systems import SingleDoFLinearSystem
        >>> sys = SingleDoFLinearSystem()
        >>> sys.constants['radius'] = 5.0
        >>> sys.constants['center_y'] = 10.0
        >>> sys.coordinates['center_x'] = 0.0
        >>> def plot(radius, center_x, center_y, time):
        ...     fig, ax = plt.subplots(1, 1)
        ...     circle = plt.Circle((center_x, center_y), radius=radius)
        ...     ax.add_patch(circle)
        ...     ax.set_title(time)
        ...     return fig, circle, ax
        ...
        >>> sys.config_plot_func = plot
        >>> sys.plot_configuration()  # doctest: +SKIP

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
        >>> from resonance.linear_systems import SingleDoFLinearSystem
        >>> sys = SingleDoFLinearSystem()
        >>> sys.constants['radius'] = 5.0
        >>> sys.constants['center_y'] = 10.0
        >>> sys.coordinates['center_x'] = 0.0
        >>> sys.speeds['center_vx'] = 0.0
        >>> def calc_coeffs():
        ...     # dummy placeholder
        ...     return 1.0, 1.0, 1.0
        >>> sys.canonical_coeffs_func = calc_coeffs
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
        ...
        >>> sys.config_plot_update_func = update
        >>> sys.free_response(1.0)  #doctest: +SKIP
              center_x  center_x_acc  center_vx
        time
        0.00       0.0           0.0        0.0
        0.01       0.0           0.0        0.0
        0.02       0.0           0.0        0.0
        0.03       0.0           0.0        0.0
        0.04       0.0           0.0        0.0
        ...        ...           ...        ...
        0.96       0.0           0.0        0.0
        0.97       0.0           0.0        0.0
        0.98       0.0           0.0        0.0
        0.99       0.0           0.0        0.0
        1.00       0.0           0.0        0.0

        [101 rows x 3 columns]
        >>> sys.animate_configuration()  #doctest: +SKIP


        """
        return self._config_plot_update_func

    @config_plot_update_func.setter
    def config_plot_update_func(self, func):
        self._config_plot_update_func = func

    def _check_system(self):

        num_coords = len(self.coordinates)
        num_speeds = len(self.speeds)

        msg = ('Check whether you have the correct number of coordinates and '
               'speeds defined. There are {} coordinates and {} speeds '
               'defined. There should be one speed for each coordinate.')

        if num_coords != num_speeds:
            raise ValueError(msg.format(num_coords, num_speeds))

    def _get_par_vals(self, par_name):
        """Returns the value of any variable stored in the parameters,
        coordinates, or measurements dictionaries."""
        # TODO : This check for duplicates is slow and called way too often.
        self._measurements._check_for_duplicate_keys()
        # TODO : This is duplicate of similar code in
        # _MeasurementsDict._compute_value(). Shouldn't have this redundancy.
        msg = '{} is not in constants, coordinates, speeds, or measurements.'
        if par_name.lower() == 'time':
            return self._time['t']
        elif par_name in self.constants:
            return self.constants[par_name]
        elif par_name in self.coordinates:
            return self.coordinates[par_name]
        elif par_name in self.speeds:
            return self.speeds[par_name]
        elif par_name in self._measurements:
            return self.measurements[par_name]
        else:
            raise KeyError(msg.format(par_name))

    def _check_meas_func(self, func):

        func_args = getfullargspec(func).args

        try:
            res = func(*np.random.random(len(func_args)))
        except Exception as e:
            msg = ("Measurement function failed to compute when given floats "
                   "as all of the arguments with: ")
            raise type(e)(msg + str(e)).with_traceback(sys.exc_info()[2])
        msg = ("Your function does not return a single float when"
               " passed floats as arguments. It returns something with type:"
               " {}. Adjust your function so that it returns a float when all"
               " arguments are floats.")
        if not isinstance(res, float):
            raise TypeError(msg.format(type(res)))

        # the function may only be a function of constants, so only do the test
        # if the measurement function has coordinates, speeds, or time.
        # TODO : Not sure if this also checks for arrays for measurements in
        # the func sig, maybe.
        has_arrays = False
        for a in func_args:
            if a in self.states or a == 'time':
                has_arrays = True

        if has_arrays:

            # constants always get a single random float and any coordinates,
            # speeds, or time get replaced by a 1D array
            arr_len = 5
            args = [random.random() if a in self.constants else
                    np.random.random(arr_len) for a in func_args]

            msg1 = ("Measurement function failed to compute when given equal "
                    "length 1D NumPy arrays as arguments for the coordinates, "
                    "speeds, measurements, and/or time and floats for the "
                    "constants with this error: ")
            msg2 = ("Your function does not return a 1D NumPy array when "
                    "passed 1D arrays of equal lengths for the coordinates, "
                    "speeds, measurements, and/or time. It returns a {}. "
                    "Adjust your function so that it returns a 1D array in "
                    "this case.")
            try:
                res = func(*args)
            except Exception as e:
                raise type(e)(msg1 + str(e)).with_traceback(sys.exc_info()[2])
            if not isinstance(res, np.ndarray):
                raise TypeError(msg2.format(type(res)))
            if res.shape != (arr_len, ):
                raise TypeError(msg2.format(type(res)))

    def add_measurement(self, name, func):
        """Creates a new measurement entry in the measurements attribute that
        uses the provided function to compute the measurement given a subset of
        the constants, coordinates, speeds, other measurements, and time.

        Parameters
        ==========
        name : string
            This must be a valid Python variable name and it should not clash
            with any names in the constants, coordinates, or speeds dictionary.
            This string can be different that the function name.
        func : function
            This function must only have existing constant, coordinate, speed,
            or measurement names, and/or the special name ``'time'`` in the
            function signature. These can be a subset of the available choices
            in constants, coordinates, speeds, measurements and any order in
            the signature is permitted. The function must be able to operate on
            both inputs that are a collection of floats or a collection of
            equal length 1D NumPy arrays and floats, i.e. the function must be
            vectorized. So be sure to use NumPy vectorized functions inside
            your function, i.e.  ``numpy.sin()`` instead of ``math.sin()``. The
            measurement function you create should return a item, either a
            scalar or array, that gives the values of the measurement.

        Examples
        ========

        >>> import numpy as np
        >>> from resonance.linear_systems import SingleDoFLinearSystem
        >>> sys = SingleDoFLinearSystem()
        >>> sys.constants['m'] = 1.0  # kg
        >>> sys.constants['c'] = 0.2  # kg*s
        >>> sys.constants['k'] = 10.0  # N/m
        >>> sys.coordinates['x'] = 1.0  # m
        >>> sys.speeds['v'] = 0.25  # m/s
        >>> def can_coeffs(m, c, k):
        ...     return m, c, k
        ...
        >>> sys.canonical_coeffs_func = can_coeffs
        >>> def force(x, v, c, k, time):
        ...    return -k * x - c * v + 5.0 * time
        ...
        >>> # The measurement function you create must be vectorized, such
        >>> # that it works with both floats and 1D arrays. For example with
        >>> # floats:
        >>> force(1.0, 0.5, 0.2, 10.0, 0.1)
        -9.6
        >>> # And with 1D arrays:
        >>> force(np.array([1.0, 1.0]), np.array([0.25, 0.25]), 0.2, 10.0,
        ...       np.array([0.1, 0.2]))
        array([-9.55, -9.05])
        >>> sys.add_measurement('f', force)
        >>> sys.measurements['f']  # time is 0.0 by default
        -10.05
        >>> sys.constants['k'] = 20.0  # N/m
        >>> sys.measurements['f']
        -20.05
        >>> # Note that you should use NumPy functions to ensure your
        >>> # measurement is vectorized.
        >>> sys.constants['force'] = 2.0
        >>> def force_mag(force):
        ...     return np.abs(force)
        ...
        >>> force_mag(-10.05)
        10.05
        >>> force_mag(np.array([-10.05, -20.05]))
        array([10.05, 20.05])
        >>> sys.add_measurement('fmag', force_mag)
        >>> sys.measurements['fmag']
        2.0

        """
        if name.lower() == 'time':
            msg = ('{} is a reserved parameter name. '
                   'Choose something different.')
            raise ValueError(msg.format(name))
        elif any([name.endswith(s) for s in _FORBIDDEN_SUFFIXES]):
            msg = ('{} are reserved suffixes. '
                   'Choose something different.')
            raise ValueError(msg.format(_FORBIDDEN_SUFFIXES))
        elif name in (list(self.constants.keys()) +
                      list(self.coordinates.keys()) +
                      list(self.speeds.keys())):
            msg = ('{} is already used as a constant or coordinate name. '
                   'Choose something different.')
            raise ValueError(msg.format(name))

        self._check_meas_func(func)

        self.measurements._funcs[name] = func
        dict.__setitem__(self.measurements, name,
                         self.measurements._compute_value(name))

    def _state_traj_to_dataframe(self, times, pos, vel, acc):

        # pads with singleton dimension if pos, vel, acc are 1D
        pos = np.atleast_2d(pos)
        vel = np.atleast_2d(vel)
        acc = np.atleast_2d(acc)

        # all should be n x m
        assert pos.shape == (len(self.coordinates), len(times))
        assert vel.shape == (len(self.speeds), len(times))
        assert acc.shape == (len(self.speeds), len(times))

        # TODO : What if they added a coordinate with the acc names?
        data = {}
        for i, c_name in enumerate(list(self.coordinates.keys())):
            data[c_name] = pos[i]
            data[c_name + self._acc_append] = acc[i]
        for i, s_name in enumerate(list(self.speeds.keys())):
            data[s_name] = vel[i]
        df = pd.DataFrame(data, index=times)
        df.index.name = self._time_var_name

        # TODO : Allow acc to be used in measurements.
        # store current values of coords, speeds, and time
        x0s = list(self.coordinates.values())
        v0s = list(self.speeds.values())
        t0 = self._time['t']

        # set coord, speeds, and time to arrays
        for i, c_name in enumerate(list(self.coordinates.keys())):
            self.coordinates[c_name] = pos[i]
        for i, s_name in enumerate(list(self.speeds.keys())):
            self.speeds[s_name] = vel[i]
        self._time['t'] = times

        # compute each measurement
        for name, value in self.measurements.items():
            df[name] = value

        # set the coords, speeds, and time back to original values
        for i, c_name in enumerate(list(self.coordinates.keys())):
            self.coordinates[c_name] = x0s[i]
        for i, s_name in enumerate(list(self.speeds.keys())):
            self.speeds[s_name] = v0s[i]
        self._time['t'] = t0

        return df

    @staticmethod
    def _calc_times(final_time, initial_time, sample_rate):
        # TODO : Should have the option to pass in unequally spaced monotonic
        # time arrays.

        if final_time < initial_time:
            raise ValueError('Final time must be greater than initial time.')

        delta = final_time - initial_time
        num = int(round(sample_rate * delta)) + 1
        return np.linspace(initial_time, final_time, num=num)

    def free_response(self, final_time, initial_time=0.0, sample_rate=100,
                      integrator="rungakutta4", **kwargs):
        """Returns a data frame with monotonic time values as the index and
        columns for each coordinate and measurement at the time value for that
        row. Note that this data frame is stored on the system as the variable
        ``.result`` until this method is called again, which will overwrite it.

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
            time, i.e. inclusive, along with times space equally based on the
            sample rate.
        integrator : string, optional
            Either ``rungakutta4`` or ``lsoda``. The ``rungakutta4`` option is
            a very simple implementation and the ``sample_rate`` directly
            affects the accuracy and quality of the result. The ``lsoda`` makes
            use of SciPy's ``odeint`` function which switches between two
            integrators for stiff and non-stiff portions of the simulation and
            is variable step so the sample rate does not affect the quality and
            accuracy of the result. This has no affect on single degree of
            freedom linear systems, as their solutions are computed
            analytically.
        **kwargs
            Keyword arguments to be passed to the underlying integration
            routine. For example, if ``integrator='lsoda'`` the keyword
            arguments for ``scipy.integrate.odeint()`` can be supplied. See
            https://docs.scipy.org/doc/scipy/reference/generated/scipy.integrate.odeint.html
            for the options.

        Returns
        =======
        df : pandas.DataFrame
            A data frame indexed by time with all of the coordinates and
            measurements as columns.

        """

        times = self._calc_times(final_time, initial_time, sample_rate)

        pos, vel, acc = self._generate_state_trajectories(times,
                                                          integrator=integrator,
                                                          **kwargs)

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
            msg = 'No plotting function has been assigned to config_plot_func.'
            raise ValueError(msg)
        else:
            args = []
            for k in getfullargspec(self.config_plot_func).args:
                if k == 'time':
                    args.append(self._time['t'])
                elif k == 'time__hist':
                    args.append(self._time['t'])
                elif k == 'time__futr':
                    args.append(self._time['t'])
                elif k.endswith('__hist'):
                    args.append(self._get_par_vals(k[:-6]))
                elif k.endswith('__futr'):
                    args.append(self._get_par_vals(k[:-6]))
                else:
                    args.append(self._get_par_vals(k))
            return self.config_plot_func(*args)

    def _resample_trajectories(self, sample_rate=60):

        trajectories = self.result.copy()

        new_times = self._calc_times(self.result.index[-1],
                                     self.result.index[0], sample_rate)
        # first and last will be same, so drop
        new_times_trunc = new_times[1:-1]

        num_cols = len(trajectories.columns)

        missing_vals = np.nan * np.ones((len(new_times_trunc), num_cols))

        nan_df = pd.DataFrame(missing_vals, columns=trajectories.columns,
                              index=new_times_trunc)

        df_with_missing = pd.concat((trajectories, nan_df)).sort_index()

        interpolated_df = df_with_missing.interpolate(method='index')

        trajectories = interpolated_df.loc[new_times]

        interval = 1000 / 60  # milliseconds

        return trajectories, interval

    def animate_configuration(self, fps=30, **kwargs):
        """Returns a matplotlib animation function based on the configuration
        plot and the configuration plot update function.

        Parameters
        ==========
        fps : integer
            The frames per second that should be displayed in the animation.
            The latest trajectory will be resampled via linear interpolation to
            create the correct number of frames. Note that the frame rate will
            depend on the CPU speed of the computer. You'll likely have to
            adjust this by trial and error to get something that matches well
            for your computer if you want the animation to run in real time.
        **kwargs
            Any extra keyword arguments will be passed to
            ``matplotlib.animation.FuncAnimation()``. The ``interval`` keyword
            argument will be ignored.

        """

        if self.config_plot_update_func is None:
            msg = ('No ploting update function has been assigned to '
                   'config_plot_update_func.')
            raise ValueError(msg)

        kwargs.pop('interval', None)  # ignore the user's supplied interval
        try:
            sample_rate = int(1.0 / np.diff(self.result.index).mean())
        except AttributeError:
            msg = ("No trajectory has been computed yet, so the animation "
                   "can't run. Run one of the response functions.")
            raise AttributeError(msg)

        fps = int(fps)
        if sample_rate != fps:
            trajectories, interval = self._resample_trajectories(fps)
        else:
            trajectories, interval = self.result, 1000 / sample_rate

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
            for k in getfullargspec(self.config_plot_update_func).args:
                if k == 'time':
                    args.append(time)
                elif k == 'time__hist':
                    args.append(trajectories[:time].index)
                elif k == 'time__futr':
                    args.append(trajectories[time:].index)
                elif k.endswith('__hist'):
                    args.append(trajectories[k[:-6]][:time])
                elif k.endswith('__futr'):
                    args.append(trajectories[k[:-6]][time:])
                elif k in trajectories:  # constants, coordinates, measurements
                    args.append(row[k])
                elif k in self.constants:
                    args.append(self.constants[k])
                else:  # must be matplotlib object
                    # TODO : This last part is quite fragile. It requires these
                    # remaining args to be in the same order as the returned
                    # tuple from the plot function and there is no way to know
                    # if these are the correct objects to append if the order
                    # isn't correct.
                    args.append(pop_list.pop(0))
            self.config_plot_update_func(*args)

        # TODO : Run this with the initial conditions so that errors will
        # propogate before the animation is run.
        # NOTE : This is useful to uncomment in debugging because the errors
        # push to the top if in the FuncAnimation.
        #gen_frame((1.0, self.result.iloc[0]), list(objs_to_modify))

        # NOTE : If the iterrows() generator is not converted to a list then
        # matplotlib will throw a StopIteration error when the animation
        # reaches the last frame instead of repeating. This causes headaches in
        # the notebook and elsewhere. See issue #39 in the resonance repo.
        return animation.FuncAnimation(fig, gen_frame,
                                       fargs=(objs_to_modify, ),
                                       frames=list(trajectories.iterrows()),
                                       interval=interval,
                                       **kwargs)
