import math
import collections
from inspect import getargspec

import numpy as np
import pandas as pd
import matplotlib.animation as animation


class ParametersDict(collections.MutableMapping, dict):

    def __getitem__(self, key):
        return dict.__getitem__(self, key)

    def __setitem__(self, key, value):
        # TODO : It would be nice to check to ensure that these names do not
        # clash with names in measurements or coordinates.
        if not key.isidentifier():
            msg = ('{} is not a valid parameter name. '
                   'It must be a valid Python variable name.')
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


class MeasurementsDict(collections.MutableMapping, dict):

    def _compute_value(self, key):

        func = self._funcs[key]

        def get_par(k):
            try:
                v = self._par[k]
            except KeyError:
                try:
                    v = self._coord[k]
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
    """
    abstract base class

    parameters : ParameterDict
    coordinates : ParameterDict
    canonical_coefficients : returns m, c, k scalars


    """

    # TODO : Only allow a single coordinate to be set on a 1 DoF system.

    def __init__(self):

        # TODO : Allow pars, coords, and meas to be set on intialization.

        self._parameters = ParametersDict({})
        self._coordinates = ParametersDict({})
        self._measurements = MeasurementsDict({})
        self._measurements._par = self._parameters
        self._measurements._coord = self._coordinates
        self._measurements._funcs = {}

        self.config_plot_func = None
        self.config_plot_update_func = None

    @property
    def parameters(self):
        return self._parameters

    @parameters.setter
    def parameters(self, val):
        raise ValueError('Not allowed to set the parameters dictionary.')

    @property
    def coordinates(self):
        return self._coordinates

    @coordinates.setter
    def coordinatess(self, val):
        raise ValueError('Not allowed to set the coordinates dictionary.')

    @property
    def measurements(self):
        return self._measurements

    @measurements.setter
    def measurementss(self, val):
        raise ValueError('Not allowed to set the measurements dictionary.')

    def _get_par_vals(self, par_name):
        """Returns the value of any variable stored in the parameters,
        coordinates, or measurements dictionaries."""
        # TODO : Maybe just merging the three dicts is better? What to do about
        # name clashes in the dictionaries? Garbage in garbage out?
        try:
            v = self.parameters[par_name]
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
            array, that gives the values of the measurement. For example::


        Examples
        ========

        >>> def f(par2, meas4, par1, coord5):
                  return par2 + meas4 + par1 + coord5
        >>> sys.add_measurement('new_meas', f)
        >>> sys.measurements['new_name']
        4.0

        """
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

    def free_response(self, initial_time, final_time, num=None):
        if num is None:
            delta = final_time - initial_time
            num = int(60 * delta)  # 60 hertz, if time is in seconds
        times = np.linspace(initial_time, final_time, num=num)

        sol_func = self._solution_func()

        coordinate_traj = sol_func(times)

        coord_name = list(self.coordinates.keys())[0]

        df = pd.DataFrame({coord_name: coordinate_traj}, index=times)
        df.index.name = 'time'

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

        args = []
        for k in getargspec(self.config_plot_func).args:
            if k == 'time':
                args.append(0.0)  # static config defaults to t=0.0
            else:
                args.append(self._get_par_vals(k))

        return self.config_plot_func(*args)

    def animate_configuration(self, **kwargs):

        fig, *objs_to_modify = self.plot_configuration()

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

    def __init__(self):

        super(BookOnCupSystem, self).__init__()

        self.parameters['height'] = 0.029  # m
        self.parameters['length'] = 0.238  # m
        self.parameters['radius'] = 0.042  # m
        self.parameters['mass'] = 1.058  # kg
        self.coordinates['book_angle'] = 0.0  # rad

    # TODO : This needs to be added in the super class with the add_coef_func()
    # method.
    def _canonical_coefficients(self):
        """A 1 DoF second order system should return the mass, damping, and
        stiffness coefficients."""

        def coeffs(height, length, radius):
            """Students will write this function themselves and pass it into
            the class via self.add_coef_func() when they get to modeling."""
            g = 9.81
            m = height**2 / 3 + length**2 / 12
            c = 0.0
            k = g * radius - g * height / 2
            return m, c, k

        args = [self._get_par_vals(k) for k in getargspec(coeffs).args]

        return coeffs(*args)
