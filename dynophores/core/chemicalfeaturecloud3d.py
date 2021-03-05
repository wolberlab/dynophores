"""
dynophores.core.chemicalfeaturecloud3d

Handles the ChemicalFeatureCloud3D class, which describes the point cloud for one superfeature.
"""


class ChemicalFeatureCloud3D:
    """
    Class to store a chemical feature 3D cloud for one superfeature.
    TODO: Add frame indices?

    Attributes
    ----------
    id : str
        Superfeature ID (= chemical feature 3D cloud ID).
    center : numpy.array
        X, Y, and Z coordiantes of geometrical center of all points in point cloud.
    points : numpy.array
        X, Y, and Z coordinates (columns) for each cloud point (rows).
    """

    def __init__(self, id, center, points):

        self.id = id
        self.center = center
        self.points = points
