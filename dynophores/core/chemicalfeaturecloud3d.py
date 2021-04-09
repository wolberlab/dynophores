"""
dynophores.core.chemicalfeaturecloud3d

Handles the ChemicalFeatureCloud3D class, which describes the point cloud for one superfeature.
"""

import pandas as pd

from dynophores.core.chemicalfeaturecloud3dpoint import ChemicalFeatureCloud3DPoint


class ChemicalFeatureCloud3D:
    """
    Class to store a chemical feature 3D cloud for one superfeature.
    TODO: Add frame indices?

    Attributes
    ----------
    center : numpy.array
        X, Y, and Z coordiantes of geometrical center of all points in point cloud.
    points : numpy.array
        X, Y, and Z coordinates (columns) for each cloud point (rows).
    """

    def __init__(self, center, points, **kwargs):

        self.center = center
        self.points = [
            point
            if isinstance(point, ChemicalFeatureCloud3DPoint)
            else ChemicalFeatureCloud3DPoint(**point)
            for point in points
        ]

    @property
    def data(self):
        """
        Superfeature cloud point data.

        Returns
        -------
        pandas.DataFrame
            Cloud point data (coordinates).
        """

        return pd.DataFrame(
            [[point.x, point.y, point.z, point.frame_ix, point.weight] for point in self.points],
            columns=["x", "y", "z", "frame_ix", "weight"],
        )
