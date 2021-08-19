"""
Unit tests for dynophores.utils.
"""

import pytest
import numpy as np

from dynophores import utils


@pytest.mark.parametrize(
    "hex, sequence_length, min_saturation, rgb_saturation_sequence",
    [
        (
            "#f73e3e",
            4,
            0.2,
            np.array(
                [
                    [0.96862745, 0.24313725, 0.24313725],
                    [0.96862745, 0.42039216, 0.42039216],
                    [0.96862745, 0.59764706, 0.59764706],
                    [0.96862745, 0.77490196, 0.77490196],
                ]
            ),
        )
    ],
)
def test_hex_to_rgb_saturation_sequence(
    hex, sequence_length, min_saturation, rgb_saturation_sequence
):
    rgb_saturation_sequence_calculated = utils.hex_to_rgb_saturation_sequence(
        hex, sequence_length, min_saturation
    )
    assert np.allclose(rgb_saturation_sequence_calculated, rgb_saturation_sequence, equal_nan=True)

@pytest.mark.parametrize(
    "hex, sequence_length, min_saturation",
    [
        ("#f73e3e", 1, 0.2),  # Sequence too short
        ("#f73e3e", 4, 10)  # Saturation not in [0, 1]
    ],
)
def test_hex_to_rgb_saturation_sequence_raise(
    hex, sequence_length, min_saturation
):
    with pytest.raises(ValueError):
        utils.hex_to_rgb_saturation_sequence(hex, sequence_length, min_saturation)