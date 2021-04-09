"""
Unit tests for dynophore.core.chemicalfeaturecloud3dpoint.ChemicalFeatureCloud3DPoint class.

Uses fixture tests.conftest.superfeature.
"""

from pathlib import Path

import numpy as np

from dynophores import parsers
from dynophores.core.chemicalfeaturecloud3d import ChemicalFeatureCloud3D
from dynophores.core.chemicalfeaturecloud3dpoint import ChemicalFeatureCloud3DPoint


PATH_TEST_DATA = Path(__name__).parent / "dynophores/tests/data"


class TestsChemicalFeatureCloud3DPoint:
    """
    Test ChemicalFeatureCloud3D class methods.
    """

    def test_init(self):

        dynophore_dict = parsers._json_pml_to_dict(
            PATH_TEST_DATA / "out/1KE7_dynophore.json",
            PATH_TEST_DATA / "out/1KE7_dynophore.pml",
        )
        superfeature_dict = next(iter(dynophore_dict["superfeatures"].values()))
        cloud_dict = superfeature_dict["cloud"]
        cloud = ChemicalFeatureCloud3D(**cloud_dict)
        point = next(iter(cloud.points))
        assert isinstance(point, ChemicalFeatureCloud3DPoint)
        assert list(point.__dict__) == ["x", "y", "z", "frame_ix", "weight"]

        # Test class attributes - check for data types
        assert isinstance(point.x, float)
        assert isinstance(point.y, float)
        assert isinstance(point.z, float)
        assert isinstance(point.frame_ix, int)
        assert isinstance(point.weight, float)

    def test_coordinates(self, chemicalfeaturecloud3dpoint):
        """
        Test class property.
        """

        assert isinstance(chemicalfeaturecloud3dpoint.coordinates, np.ndarray)
