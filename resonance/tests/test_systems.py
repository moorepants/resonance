import pytest

from ..system import (_ConstantsDict, _MeasurementsDict, _CoordinatesDict,
                      _StatesDict)


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
