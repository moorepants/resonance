import pytest

from ..linear_systems import *


def test_nonvalid_parameters_key():
    p = ParametersDict({})
    with pytest.raises(ValueError):
        p['ben stiller'] = 12.0


def test_setting_measurements_item():
    m = MeasurementsDict({})
    with pytest.raises(ValueError):
        m['a'] = 12.0
