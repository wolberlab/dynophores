"""
Unit tests for dynophore.core.envpartner.EnvPartner class.

Uses fixture tests.conftest.envpartner.
"""

from pathlib import Path

import pytest
import numpy as np

from dynophores import parsers
from dynophores.core.envpartner import EnvPartner

PATH_TEST_DATA = Path(__name__).parent / "dynophores/tests/data"


class TestsEnvPartner:
    """
    Test EnvPartner class methods.
    """

    def test_init(self):

        dynophore_dict = parsers._json_pml_to_dict(
            PATH_TEST_DATA / "out/1KE7_dynophore.json",
            PATH_TEST_DATA / "out/1KE7_dynophore.pml",
        )
        superfeature_dict = next(iter(dynophore_dict["superfeatures"].values()))
        envpartner_dict = next(iter(superfeature_dict["envpartners"].values()))
        envpartner = EnvPartner(**envpartner_dict)
        assert isinstance(envpartner, EnvPartner)
        assert list(envpartner.__dict__) == [
            "id",
            "residue_name",
            "residue_number",
            "chain",
            "atom_numbers",
            "occurrences",
            "distances",
        ]

        # Test class attributes - check for data types
        assert isinstance(envpartner.id, str)
        assert isinstance(envpartner.residue_name, str)
        assert isinstance(envpartner.residue_number, int)
        assert isinstance(envpartner.chain, str)
        assert isinstance(envpartner.atom_numbers[0], int)
        assert isinstance(envpartner.occurrences[0], int)
        assert isinstance(envpartner.distances[0], float)

    @pytest.mark.parametrize(
        "envpartner_dict",
        [
            {
                "id": "ILE-10-A[169,171,172]",
                "residue_name": "ILE",
                "residue_number": 10,
                "chain": "A",
                "atom_numbers": [169, 171, 172],
                "occurrences": np.array([0, 0, 1, 1]),  # Length occurrences != distances
                "distances": np.array([6.0, 6.0, 3.0]),
            }
        ],
    )
    def test_init_raises(self, envpartner_dict):

        with pytest.raises(ValueError):
            EnvPartner(**envpartner_dict)

    @pytest.mark.parametrize("residue_id", ["ILE-10-A"])
    def test_residue_id(self, envpartner, residue_id):
        """
        Test class property.
        """

        assert envpartner.residue_id == residue_id

    @pytest.mark.parametrize("n_frames", [1002])
    def test_n_frames(self, envpartner, n_frames):
        """
        Test class property.
        """

        assert envpartner.n_frames == n_frames

    @pytest.mark.parametrize("count", [995])
    def test_count(self, envpartner, count):
        """
        Test class property.
        """

        assert envpartner.count == count

    @pytest.mark.parametrize("frequency", [99.3])
    def test_frequency(self, envpartner, frequency):
        """
        Test class property.
        """

        assert envpartner.frequency == frequency
