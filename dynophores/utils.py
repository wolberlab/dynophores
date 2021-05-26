"""
Utility functions.
"""


def hex_to_rgb(hex_string, scale=True):
    """
    Convert hex code to RGB code.

    Parameters
    ----------
    hex_string : str
        6-letter hex string.

    Return
    ------
    list of int or float
        RGB code, e.g. (255, 255, 255).
    """

    if not isinstance(hex_string, str):
        raise ValueError(f"Input must be string but is {type(hex_string)}.")
    hex_string = hex_string.lstrip("#")
    if len(hex_string) != 6:
        raise ValueError(f"Input hex string must be of length 6 but is '{hex_string}'.")
    r_hex = hex_string[0:2]
    g_hex = hex_string[2:4]
    b_hex = hex_string[4:6]

    rgb = [int(r_hex, 16), int(g_hex, 16), int(b_hex, 16)]

    if scale:
        rgb = [round(value / float(255), 3) for value in rgb]

    return rgb
