"""
Unit tests for dynophore.core.superfeature.SuperFeature class.

Uses fixture tests.conftest.superfeature.
"""

from pathlib import Path

import pytest
import pandas as pd

from dynophores import parsers
from dynophores.core.superfeature import SuperFeature
from dynophores.core.envpartner import EnvPartner
from dynophores.core.chemicalfeaturecloud3d import ChemicalFeatureCloud3D


PATH_TEST_DATA = Path(__name__).parent / "dynophores/tests/data"


class TestsSuperFeature:
    """
    Test SuperFeature class methods.
    """

    def test_init(self):

        dynophore_dict = parsers._json_pml_to_dict(
            PATH_TEST_DATA / "out/1KE7_dynophore.json",
            PATH_TEST_DATA / "out/1KE7_dynophore.pml",
        )
        superfeature_dict = next(iter(dynophore_dict["superfeatures"].values()))
        superfeature = SuperFeature(**superfeature_dict)
        assert isinstance(superfeature, SuperFeature)
        assert list(superfeature.__dict__) == [
            "id",
            "feature_type",
            "atom_numbers",
            "occurrences",
            "envpartners",
            "color",
            "cloud",
        ]

        # Test class attributes - check for data types
        assert isinstance(superfeature.id, str)
        assert isinstance(superfeature.feature_type, str)
        assert isinstance(superfeature.atom_numbers[0], int)
        assert isinstance(superfeature.occurrences[0], int)
        assert isinstance(superfeature.envpartners, dict)
        assert isinstance(next(iter(superfeature.envpartners.values())), EnvPartner)
        assert isinstance(superfeature.color, str)
        assert isinstance(superfeature.cloud, ChemicalFeatureCloud3D)

    def test_envpartners_occurrences(self, superfeature):
        """
        Test class property.
        """

        data = superfeature.envpartners_occurrences
        assert isinstance(data, pd.DataFrame)
        assert data.index.to_list() == list(range(0, len(superfeature.occurrences)))
        assert sorted(data.columns.to_list()) == sorted(
            [envpartner.id for _, envpartner in superfeature.envpartners.items()]
        )
        assert data.dtypes.unique() == "int32"

    
    @pytest.mark.parametrize("envpartners_collapsed", [
        ["ILE-10-A[169,171,172]", "PHE-82-A[1245,1246,1247,1248,1249,1250]"]
        ]
    )
    def test_envpartners_occurrences_collapsed(self, superfeature, envpartners_collapsed):
        """
        Test class property.
        """

        data = superfeature.envpartners_occurrences_collapsed
        assert isinstance(data, pd.DataFrame)
        assert data.index.to_list() == list(range(0, len(superfeature.occurrences)))
        assert sorted(data.columns.to_list()) == sorted(envpartners_collapsed)
        assert data.dtypes.unique() == "int32"

    def test_envpartners_distances(self, superfeature):
        """
        Test class property.
        """

        data = superfeature.envpartners_distances
        assert isinstance(data, pd.DataFrame)
        assert data.index.to_list() == list(range(0, len(superfeature.occurrences)))
        assert sorted(data.columns.to_list()) == sorted(
            [envpartner.id for _, envpartner in superfeature.envpartners.items()]
        )
        assert data.dtypes.unique() == "float64"

    @pytest.mark.parametrize("data_type", ["xxx"])
    def test_data_raises(self, superfeature, data_type):
        """
        Test helper method to call envpartners_occurrences and envpartners_distances properties.
        """

        with pytest.raises(KeyError):
            superfeature._data(data_type)

    @pytest.mark.parametrize("n_frames", [1002])
    def test_n_frames(self, superfeature, n_frames):
        """
        Test class property.
        """

        assert superfeature.n_frames == n_frames

    @pytest.mark.parametrize(
        "count, frequency, envpartner_ids",
        [
            (
                [1001, 995, 945, 57],
                [99.90, 99.30, 94.31, 5.69],
                [
                    "any",
                    "ILE-10-A[169,171,172]",
                    "ILE-10-A[169,171]",
                    "PHE-82-A[1245,1246,1247,1248,1249,1250]",
                ],
            )
        ],
    )
    def test_count_frequency(self, superfeature, count, frequency, envpartner_ids):
        """
        Test class property.
        """

        count = pd.Series(count, index=envpartner_ids)
        assert all(superfeature.count == count)

        frequency = pd.Series(frequency, index=envpartner_ids)
        assert all(superfeature.frequency == frequency)
