"""
dynophores.tests.fixures

Fixures to be used in unit testing.
"""

from pathlib import Path

import pytest
import numpy as np

from dynophores import Dynophore, SuperFeature, EnvPartner, ChemicalFeatureCloud3D

PATH_TEST_DATA = Path(__name__).parent / "dynophores" / "tests" / "data"


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
def chemicalfeaturecloud3d():

    cloud_dict = {
        "id": "H[4599,4602,4601,4608,4609,4600]",
        "center": np.array([1, 1, 1]),
        "points": np.array([[-1, -1, -1], [0, 0, 0], [1, 1, 1]]),
    }
    cloud = ChemicalFeatureCloud3D(**cloud_dict)
    return cloud


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

    cloud_dict = {
        "id": "H[4599,4602,4601,4608,4609,4600]",
        "center": np.array([1, 1, 1]),
        "points": np.array([[-1, -1, -1], [0, 0, 0], [1, 1, 1]]),
    }

    superfeature_dict = {
        "id": "H[4599,4602,4601,4608,4609,4600]",
        "feature_type": "H",
        "atom_numbers": [4599, 4602, 4601, 4608, 4609, 4600],
        "occurrences": np.array([0, 1, 1]),
        "envpartners": {
            envpartner1_dict["id"]: EnvPartner(**envpartner1_dict),
            envpartner2_dict["id"]: EnvPartner(**envpartner2_dict),
        },
        "cloud": ChemicalFeatureCloud3D(**cloud_dict),
    }

    superfeature = SuperFeature(**superfeature_dict)
    return superfeature


@pytest.fixture(scope="module")
def dynophore():

    dynophore = Dynophore.from_dir(PATH_TEST_DATA / "out")
    return dynophore
