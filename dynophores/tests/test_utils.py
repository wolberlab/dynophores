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
            3,
            0,
            np.array(
                [[0, 0, 0], [0, 0, 0], [0, 0, 0]],
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
    assert np.array_equal(rgb_saturation_sequence_calculated, rgb_saturation_sequence)
