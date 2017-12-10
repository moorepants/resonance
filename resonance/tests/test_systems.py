from math import isclose

import numpy as np
import matplotlib.pyplot as plt
import pytest

from ..system import (System, _ConstantsDict, _MeasurementsDict,
                      _SingleDoFCoordinatesDict, _StatesDict)


def test_nonvalid_parameters_key():
    p = _ConstantsDict({})
    with pytest.raises(ValueError):
        p['ben stiller'] = 12.0
    with pytest.raises(ValueError):
        p['time'] = 12.0
    with pytest.raises(ValueError):
        p['time__hist'] = 12.0
    p['a'] = 1.0
    if 'a' in p:
        pass
    del p['a']


def test_setting_measurements_item():
    m = _MeasurementsDict({})
    # not allowed to set measurements
    with pytest.raises(ValueError):
        m['a'] = 12.0


def test_setting_coordinates_item():
    c = _SingleDoFCoordinatesDict({})

    # not a python identifier
    with pytest.raises(ValueError):
        c['a '] = 12.0

    # reserved variable name
    with pytest.raises(ValueError):
        c['time'] = 12.0

    # reserved variable name
    with pytest.raises(ValueError):
        c['time__hist'] = 12.0

    # can only add one coordinate
    c['first_key'] = 12.0
    with pytest.raises(ValueError):
        c['second_key'] = 12.0

    assert 'first_key' in c

    del c['first_key']

    assert 'first_key' not in c


def test_setting_state_item():
    s = _StatesDict({})
    with pytest.raises(ValueError):
        s['a'] = 0.0


def test_system():
    sys = System()
    sys.constants['a'] = 1.0
    assert isclose(sys._time['t'], 0.0)
    assert isclose(sys._get_par_vals('time'), 0.0)

    def plot(a):
        fig, ax = plt.subplots(1, 1)
        return fig

    def update(a):
        pass

    sys.config_plot_func = None
    sys.config_plot_update_func = None
    with pytest.raises(ValueError):
        sys.animate_configuration()

    sys.config_plot_func = plot
    sys.config_plot_update_func = None
    with pytest.raises(ValueError):
        sys.animate_configuration()

    # no trajectory computed
    sys.config_plot_func = plot
    sys.config_plot_update_func = update
    with pytest.raises(AttributeError):
        sys.animate_configuration()

    # no speed is set
    sys.coordinates['x'] = 1.0
    with pytest.raises(ValueError):
        sys._check_system()


def test_add_measurement():
    sys = System()
    sys.constants['a'] = 1.0
    sys.constants['b'] = 2.0
    sys.coordinates['x'] = 3.0
    sys.speeds['v'] = 4.0

    # this will fail if all floats are passed in
    def meas1(a, b, x, v, time):
        return a * np.array([x + v])

    with pytest.raises(TypeError):
        sys.add_measurement('meas1', meas1)

    # this will fail if arrays are passed in for x, v, time
    def meas2(a, b, x, v, time):
        return np.sum([x, v, time])

    with pytest.raises(TypeError):
        sys.add_measurement('meas2', meas2)

    # this will fail because the function just doesn't work
    def meas3(a, b, x, v, time):
        return str(a) + int(b)

    with pytest.raises(TypeError):
        sys.add_measurement('meas3', meas3)

    def meas4(a, b, x, v, time):
        return a**2 * v * time

    sys.add_measurement('meas4', meas4)

    def meas5(a, b, x, v, time, meas4):
        return np.sum(meas4)

    with pytest.raises(TypeError):
        sys.add_measurement('meas5', meas5)

    del sys.measurements['meas4']

    # should be able to add only a function of constants too (and array check
    # doesn't happen)
    def meas6(a, b):
        return a + b

    sys.add_measurement('meas6', meas6)
