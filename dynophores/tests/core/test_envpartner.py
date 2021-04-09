"""
Unit tests for dynophore.core.envpartner.EnvPartner class.

Uses fixture tests.conftest.envpartner.
"""

import pytest
import numpy as np

from dynophores.core.envpartner import EnvPartner


class TestsEnvPartner:
    """
    Test EnvPartner class methods.
    """

    @pytest.mark.parametrize(
        "envpartner_dict",
        [
            {
                "id": "ILE-10-A[169,171,172]",
                "residue_name": "ILE",
                "residue_number": 10,
                "chain": "A",
                "atom_numbers": [169, 171, 172],
                "occurrences": np.array([0, 0, 1]),
                "distances": np.array([6.0, 6.0, 3.0]),
            }
        ],
    )
    def test_init(self, envpartner_dict):

        envpartner = EnvPartner(**envpartner_dict)
        assert isinstance(envpartner, EnvPartner)

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

    @pytest.mark.parametrize("n_frames", [3])
    def test_n_frames(self, envpartner, n_frames):

        assert envpartner.n_frames == n_frames

    @pytest.mark.parametrize("count", [1])
    def test_count(self, envpartner, count):

        assert envpartner.count == count

    @pytest.mark.parametrize("frequency", [1 / 3.0])
    def test_frequency(self, envpartner, frequency):

        frequency = round(frequency * 100, 2)
        assert envpartner.frequency == frequency
