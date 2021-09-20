"""
Handles the ChemicalFeatureCloud3DPoint class, which describes a point in a superfeature cloud.
"""

import numpy as np


class ChemicalFeatureCloud3DPoint:
    """
    Class to store a point in a superfeature cloud.

    Attributes
    ----------
    x : float
        X coordinate of the point.
    y : float
        Y coordinate of the point.
    z : float
        Z coordinate of the point.
    frame_ix : int
        Frame index in which the cloud point occurs.
    weight : float
        Cloud point weight.
    """

    def __init__(self, x, y, z, frame_ix, weight, **kwargs):

        self.x = x
        self.y = y
        self.z = z
        self.frame_ix = frame_ix
        self.weight = weight

    @property
    def coordinates(self):
        return np.array([self.x, self.y, self.z])
