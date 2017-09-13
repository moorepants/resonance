import collections
from inspect import getargspec

import numpy as np
import pandas as pd
import matplotlib.animation as animation


class MeasurmentsDict(collections.MutableMapping, dict):

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
        print('Warning : cannot set measurement values.')
        # ignore and override the value
        dict.__setitem__(self, key, self._compute_value(key))

    def __delitem__(self, key):
        dict.__delitem__(self, key)

    def __iter__(self):
        return dict.__iter__(self)

    def __len__(self):
        return dict.__len__(self)

    def __contains__(self, x):
        return dict.__contains__(self, x)


class BookOnCupSystem(object):
    """This will eventually be a subclass of a 1 DoF canoncial form second
    order linear system."""
    def __init__(self):
        self.parameters = {'height': 0.029,  # m
                           'length': 0.238,  # m
                           'radius': 0.042,  # m
                           'mass': 1.058}  # kg

        self.coordinates = {'book_angle': 0.0}

        # include references to the necessary attributes on the measurement
        # dict, or should I just store a reference to self on the measurement
        # dict?
        self.measurements = MeasurmentsDict({})
        self.measurements._funcs = {}
        self.measurements._par = self.parameters
        self.measurements._coord = self.coordinates

        self.configuration_plot_function = None
        self.configuration_plot_update_function = None

    def _get_par_vals(self, par_name):
        # TODO : Maybe just merging the three dicts is better?
        try:
            v = self.parameters[par_name]
        except KeyError:
            try:
                v = self.coordinates[par_name]
            except KeyError:
                v = self.measurements[par_name]
        return v

    # TODO : Should this be moved to the setitem method of the MeasurementDict?
    # Then you could do: sys.measurements['distance'] = lambda x1, x2: x1 + x2.
    # I'm not sure if that would be more confusing, because you would set it to
    # a function and it would return a value which is unexpected.
    def add_measurement(self, name, func):
        self.measurements._funcs[name] = func
        self.measurements[name] = 0.0  # will get converted to correct number

    def _canonical_coefficients(self):
        """A 1 DoF second order system should return the mass, damping, and
        stiffness coefficients."""

        def coeffs(height, length, radius):
            """Students will write this function themselves and pass it into
            the class when they get to modeling."""
            g = 9.81
            m = height**2 / 3 + length**2 / 12
            c = 0.0
            k = g * radius - g * height / 2
            return m, c, k

        args = [self._get_par_vals(k) for k in getargspec(coeffs).args]

        return coeffs(*args)

    @staticmethod
    def _natural_frequency(mass, stiffness):
        return np.sqrt(stiffness / mass)

    @staticmethod
    def _damping_ratio(mass, damping, natural_frequency):
        return damping / 2.0 / mass / natural_frequency

    @staticmethod
    def _damped_natural_frequency(natural_frequency, damping_ratio):
        return natural_frequency * np.sqrt(1 - damping_ratio**2)

    def _solution_func(self):

        m, c, k = self._canonical_coefficients()

        if c == 0.0:
            sol_func = self._no_damping_solution
        else:  # damping, so check zeta
            omega_n = self._natural_frequency(m, k)
            zeta = self._damping_ratio(m, c, omega_n)
            if zeta < 1.0:
                sol_func = self._underdamped_solution
            elif zeta > 1.0:
                sol_func = self._overdamped_solution
            elif zeta == 0.0:
                sol_func = self._critically_damped_solution
            else:
                sol_func = None

        return sol_func

    def _no_damping_solution(self, time):

        print('Simulating system with no damping.')

        t = time

        m, c, k = self._canonical_coefficients()

        wn = self._natural_frequency(m, k)

        x0 = self.coordinates['book_angle']
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

        x0 = self.coordinates['book_angle']
        v0 = 0.0

        A = np.sqrt(((v0 + z * wn * x0)**2 + (x0 * wd)**2) / wd**2)
        phi = np.atan(x0 * wd / (v0 + z * wn * x0))

        return A * np.exp(-z * wn * t) * np.sin(wd * t + phi)

    def _overdamped_solution(self, time):

        t = time

        m, c, k = self._canonical_coefficients()

        wn = self._natural_frequency(m, k)
        z = self._damping_ratio(m, c, wn)

        x0 = self.coordinates['book_angle']
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

        x0 = self.coordinates['book_angle']
        v0 = 0.0

        a1 = x0
        a2 = v0 + wn * x0

        return (a1 + a2 * t) * np.exp(-wn * t)

    def free_response(self, initial_time, final_time, num=None):
        if num is None:
            delta = final_time - initial_time
            num = int(60 * delta)
        times = np.linspace(initial_time, final_time, num=num)

        sol_func = self._solution_func()

        coordinate_traj = sol_func(times)
        # TODO : Make functions that compute the derivative of the coordinate
        # solution and add the velocity trajectory.

        df = pd.DataFrame({'book_angle': coordinate_traj}, index=times)
        df.index.name = 'Time [s]'

        # TODO : Need a way to compute the measurements as array values based
        # on the book_angle changing at each time but the current method of
        # letting the measurement be computed by the stored coordinate scalar
        # is a bit problematic.
        for k, v in self.measurements.items():
            vals = np.zeros_like(times)
            x0 = self.coordinates['book_angle']
            for i, thetai in enumerate(coordinate_traj):
                self.coordinates['book_angle'] = thetai
                vals[i] = self.measurements[k]
            self.coordinates['book_angle'] = x0
            df[k] = vals

        self.result = df

        return df

    def plot_configuration(self):

        args = []
        for k in getargspec(self.configuration_plot_function).args:
            if k == 'time':
                args.append(0.0)
            else:
                args.append(self._get_par_vals(k))

        return self.configuration_plot_function(*args)

    def animate_configuration(self, **kwargs):

        fig, *objs_to_modify = self.plot_configuration()

        def gen_frame(row_tuple, pop_list):
            time = row_tuple[0]
            row = row_tuple[1]
            # Don't mutate the orginal list.
            pop_list = pop_list.copy()
            args = []
            for k in getargspec(self.configuration_plot_update_function).args:
                print(k)
                if k == 'time':
                    args.append(time)
                else:
                    try:
                        args.append(row[k])
                    except KeyError:
                        # requires these to be in the same order
                        args.append(pop_list.pop(0))
            self.configuration_plot_update_function(*args)

        return animation.FuncAnimation(fig, gen_frame, fargs=(objs_to_modify, ),
                                       frames=self.result.iterrows(), **kwargs)
