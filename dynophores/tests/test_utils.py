"""
Unit tests for dynophores.utils.
"""

import pytest

from dynophores import utils


@pytest.mark.parametrize(
    "hex_string, scale, rgb",
    [
        ("ffc20e", True, [1.0, 0.761, 0.055]),
        ("ffc20e", False, [255, 194, 14]),
        ("#ffc20e", False, [255, 194, 14]),
    ],
)
def test_hex_to_rgb(hex_string, scale, rgb):
    """
    Test hex to RGB conversion.
    """

    assert utils.hex_to_rgb(hex_string, scale) == rgb


@pytest.mark.parametrize(
    "hex_string",
    ["a", 0],
)
def test_hex_to_rgb_raises(hex_string):
    """
    Test hex to RGB conversion.
    """

    with pytest.raises(ValueError):
        assert utils.hex_to_rgb(hex_string)
