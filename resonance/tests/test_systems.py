from math import isclose

import matplotlib.pyplot as plt
import pytest

from ..system import (System, _ConstantsDict, _MeasurementsDict,
                      _CoordinatesDict, _StatesDict)


def test_nonvalid_parameters_key():
    p = _ConstantsDict({})
    with pytest.raises(ValueError):
        p['ben stiller'] = 12.0
    with pytest.raises(ValueError):
        p['time'] = 12.0
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
    m = _CoordinatesDict({})
    with pytest.raises(ValueError):
        m['a '] = 12.0

    with pytest.raises(ValueError):
        m['time'] = 12.0

    m['first_key'] = 12.0
    with pytest.raises(ValueError):
        m['second_key'] = 12.0

    if 'first_key' in m:
        pass
    del m['first_key']


def test_setting_state_item():
    s = _StatesDict({})
    with pytest.raises(ValueError):
        s['a'] = 0.0


def test_system():
    sys = System()
    sys.constants['a'] = 1.0
    assert isclose(sys._time, 0.0)
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
