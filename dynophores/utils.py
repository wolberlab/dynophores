"""
Contains utility functions.
"""

import numpy as np
from matplotlib import colors


def hex_to_rgb_saturation_sequence(hex, sequence_length, min_saturation=0.2):
    """
    Generate a sequence of RGB values from an input hex color that changes in saturation from high
    saturation at the beginning (100%) to low saturation at the end (by default 20%).

    Parameters
    ----------
    hex : str
        Hex string for input color.
    sequence_length : int
        Lenght of output sequence.
    min_saturation : float
        Minimum saturation value that can be used to avoid running into grey (default: 0.2).
    """

    if not 0 <= min_saturation <= 1:
        raise ValueError(f"Saturation value must be in [0, 1] but is {min_saturation}.")

    # Hex to (scaled) RGB
    color_rgb = colors.hex2color(hex)
    color_rgb = np.array(color_rgb)
    # RGB to HSV
    color_hsv = colors.rgb_to_hsv(color_rgb)
    # HSV sequence
    hsv_saturation_sequence = np.array(
        [
            [color_hsv[0], (1 - i / sequence_length) * min_saturation, color_hsv[2]]
            for i in range(1, sequence_length + 1)
        ]
    )
    # HSV sequence to RGB sequence
    rgb_saturation_sequence = np.apply_along_axis(colors.hsv_to_rgb, 1, hsv_saturation_sequence)
    rgb_saturation_sequence = rgb_saturation_sequence
    return rgb_saturation_sequence
