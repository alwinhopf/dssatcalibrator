"""Optional plugins must not silently fall back to DSSAT."""
import pytest
from dssatcalibrator.adapters import model_adapter, model_platform


def test_dssat_requires_no_extension():
    assert model_adapter({}) is None
    assert model_adapter({'model': {'platform': 'dssat_csm'}}) is None


def test_unknown_adapter_fails_closed():
    with pytest.raises(ValueError, match='installed adapter'):
        model_adapter({'model': {'platform': 'not-a-crop-model'}})


@pytest.mark.parametrize('name', ['apsim', 'apsim-ng', 'apsimx', 'apsim_next_generation'])
def test_platform_aliases(name):
    assert model_platform({'model': {'platform': name}}) == 'apsim_ng'
