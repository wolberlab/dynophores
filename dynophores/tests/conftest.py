"""
dynophores.tests.fixures

Fixures to be used in unit testing.
"""

import pytest
import numpy as np

from dynophores import SuperFeature, EnvPartner


@pytest.fixture(scope="module")
def envpartner():

    envpartner_dict = {
        "id": "ILE-10-A[169,171,172]",
        "residue_name": "ILE",
        "residue_number": 10,
        "chain": "A",
        "atom_numbers": [169, 171, 172],
        "occurrences": np.array([0, 0, 1]),
        "distances": np.array([6.0, 6.0, 3.0]),
    }

    envpartner = EnvPartner(**envpartner_dict)
    return envpartner


@pytest.fixture(scope="module")
def superfeature():

    envpartner1_dict = {
        "id": "ILE-10-A[169,171,172]",
        "residue_name": "ILE",
        "residue_number": 10,
        "chain": "A",
        "atom_numbers": [169, 171, 172],
        "occurrences": np.array([0, 0, 1]),
        "distances": np.array([6.0, 6.0, 3.0]),
    }

    envpartner2_dict = {
        "id": "VAL-18-A[275,276,277]",
        "residue_name": "VAL",
        "residue_number": 18,
        "chain": "A",
        "atom_numbers": [275, 276, 277],
        "occurrences": np.array([0, 1, 0]),
        "distances": np.array([6.0, 4.0, 4.0]),
    }

    superfeature_dict = {
        "id": "H[4599,4602,4601,4608,4609,4600]",
        "feature_type": "H",
        "atom_numbers": [4599, 4602, 4601, 4608, 4609, 4600],
        "occurrences": np.array([0, 1, 1]),
        "envpartners": [EnvPartner(**envpartner1_dict), EnvPartner(**envpartner2_dict)],
    }

    superfeature = SuperFeature(**superfeature_dict)
    return superfeature
