"""
Unit tests for dynophore.core.chemicalfeaturecloud3d.ChemicalFeatureCloud3D class.

Uses fixture tests.conftest.superfeature.
"""

from pathlib import Path

import numpy as np

from dynophores import parsers
from dynophores.core.chemicalfeaturecloud3d import ChemicalFeatureCloud3D
from dynophores.core.chemicalfeaturecloud3dpoint import ChemicalFeatureCloud3DPoint


PATH_TEST_DATA = Path(__name__).parent / "dynophores/tests/data"


class TestsChemicalFeatureCloud3D:
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
        assert isinstance(cloud, ChemicalFeatureCloud3D)
        assert list(cloud.__dict__) == ["center", "points"]

        # Test class attributes - check for data types
        assert isinstance(cloud.center, np.ndarray)
        assert isinstance(cloud.points[0], ChemicalFeatureCloud3DPoint)

    def test_data(self, chemicalfeaturecloud3d):
        """
        Test class property.
        """

        assert chemicalfeaturecloud3d.data.columns.to_list() == [
            "x",
            "y",
            "z",
            "frame_ix",
            "weight",
        ]
