"""
Contains utility functions.
"""

import logging
from pathlib import Path
import urllib

import numpy as np
from matplotlib import colors

_logger = logging.getLogger(__name__)


def _download_file(target_url, target_path):
    """
    Download a file from target URL to be saved as defined in the target path.

    Parameters
    ----------
    target_url : str
        Target URL.
    target_path : str or pathlib.Path
        Path to target file.
    """

    target_path = Path(target_path)

    if target_path.name != Path(target_url).name:
        raise ValueError("File name in target URL and target path must be the same.")

    if Path(target_path).exists():
        _logger.debug("Dynophore installer already available. Continue.")
    else:
        _logger.debug(f"Download from: {target_url}")
        _logger.debug(f"Download to: {target_path.parent}")
        urllib.request.urlretrieve(target_url, target_path)


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

    if sequence_length < 2:
        raise ValueError(f"Sequence must have at least 2 elements but has {sequence_length}.")
    if not 0 <= min_saturation <= 1:
        raise ValueError(f"Saturation value must be in [0, 1] but is {min_saturation}.")

    # Hex to (scaled) RGB
    color_rgb = colors.hex2color(hex)
    color_rgb = np.array(color_rgb)

    # RGB to HSV [hue, saturation, value]
    color_hsv = colors.rgb_to_hsv(color_rgb)

    # HSV sequence
    # We want to generate a sequence of colors in HSV format with descending saturation
    # (colors shall grow paler):
    # H: Keep fixed
    # S: Decrease (XXX)
    # V: Keep fixed
    hsv_saturation_sequence = []
    # Add first colors = input color
    hsv_saturation_sequence.append(color_hsv)
    # Define saturation width: Maximum (input HSV saturation) to minimum (default 0.2)
    s_width = color_hsv[1] - min_saturation
    # Define saturation step size
    s_step_size = s_width / (sequence_length - 1)
    # Add decreasing colors
    for i in range(1, sequence_length):
        h = color_hsv[0]
        # Define step size (maximum: , minimum: )
        s = color_hsv[1] - i * s_step_size
        v = color_hsv[2]
        hsv_saturation_sequence.append([h, s, v])
    hsv_saturation_sequence = np.array(hsv_saturation_sequence)

    # HSV sequence to RGB sequence
    rgb_saturation_sequence = np.apply_along_axis(colors.hsv_to_rgb, 1, hsv_saturation_sequence)
    rgb_saturation_sequence = rgb_saturation_sequence
    return rgb_saturation_sequence
