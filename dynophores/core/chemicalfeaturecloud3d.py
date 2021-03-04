"""
dynophores.core.chemicalfeaturecloud3d

Handles the ChemicalFeatureCloud3D class, which describes the point cloud for one superfeature.
"""

import pandas as pd


class ChemicalFeatureCloud3D:
    """
    Class to store a chemical feature 3D cloud for one superfeature.
    TODO: Add frame indices?

    Attributes
    ----------
    id : str
        Superfeature ID (= chemical feature 3D cloud ID).
    center : numpy.array
        Cloud center.
    coordinates : numpy.array
        X, Y, and Z coordinates (columns) for each cloud point (rows).
    """

    def __init__(self, id, center, coordinates):

        self.id = id
        self.center = center
        self._coordinates = coordinates

    @property
    def coordinates(self):
        return pd.DataFrame(self._coordinates, columns=["x", "y", "z"])
