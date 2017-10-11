from inspect import getargspec

import numpy as np
import scipy as sp
import scipy.integrate  # scipy doesn't import automatically

from .system import System as _System


class NonLinearSystem(_System):

    def _integrate_equations_of_motion(self, times, integrator='lsoda'):

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

    def _generate_state_trajectories(self, times):
        """This method should return arrays for position, velocity, and
        acceleration of the coordinates."""
        int_res = self._integrate_equations_of_motion(times)
        f = self.equations_of_motion
        args = tuple([self._get_par_vals(k) for k in getargspec(f).args[2:]])
        res = np.asarray(f(int_res.T, times, *args))
        return int_res[:, 0], int_res[:, 1], res[1, :]


class SimplePendulum(NonLinearSystem):
    pass
