"""
test the gvec.config module and related classes.
"""

import pytest

try:
    import gvec
except ImportError:
    pass  # tests will be skipped via the `check_import` fixture


# === Fixtures === #


@pytest.fixture(params=["picard_current"])
def default_config(request):
    if request.param == "picard_current":
        return gvec.core.minimize.PicardCurrentConfig()


# === Tests === #


def test_PicardCurrentConfig_defaults():
    """Test that the PicardCurrentConfig class has the correct default values."""
    config = gvec.core.minimize.PicardCurrentConfig()
    assert config.max_updates == 50
    assert config.iota_tol == 1e-8


def test_defaults(default_config):
    """Test that the default values of the config classes are correct."""
    assert len(set(default_config.keys())) > 0
    for key in default_config.keys():
        assert default_config[key] == default_config.default(key)
        assert getattr(default_config, key) == default_config[key]
        assert default_config.to_dict()[key] == default_config[key]
    assert set(default_config.to_dict()) == set(default_config.keys())
    assert type(default_config.to_dict()) is dict
