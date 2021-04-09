"""
Unit tests for dynophore.viz.view3d.parser.
"""

from collections import OrderedDict
from pathlib import Path

import pytest
import numpy as np

from dynophores import parsers

PATH_TEST_DATA = Path(__name__).parent / "dynophores" / "tests" / "data"


@pytest.mark.parametrize(
    "filepath, cloud_keys",
    [
        (
            PATH_TEST_DATA / "out/1KE7_dynophore.pml",
            ["id", "color", "center", "points"],
        )
    ],
)
def test_pml_to_dict(filepath, cloud_keys):
    """
    Test data types (not data values) in dict generated from PML file.
    """

    pml_dict = parsers._pml_to_dict(filepath)
    assert isinstance(pml_dict, OrderedDict)

    for superfeature_id, data in pml_dict.items():
        assert isinstance(data, OrderedDict)
        assert list(data.keys()) == cloud_keys
        assert isinstance(superfeature_id, str)
        assert data["id"] == superfeature_id
        assert isinstance(data["color"], str)
        assert isinstance(data["center"], np.ndarray)
        assert isinstance(data["points"], np.ndarray)


@pytest.mark.parametrize(
    "filepath, dynophore_keys, superfeature_keys, envpartner_keys",
    [
        (
            PATH_TEST_DATA / "out/1KE7_dynophore.json",
            ["id", "ligand_name", "ligand_smiles", "superfeatures"],
            ["id", "feature_type", "atom_numbers", "occurrences", "envpartners"],
            ["id", "name", "atom_numbers", "occurrences", "distances"],
        )
    ],
)
def test_json_to_dict(
    filepath,
    dynophore_keys,
    superfeature_keys,
    envpartner_keys,
):
    """
    Test data types (not data values) in dict generated from JSON file.
    """

    json_dict = parsers._json_to_dict(filepath)
    assert list(json_dict.keys()) == dynophore_keys
    assert isinstance(json_dict["superfeatures"], list)

    for superfeature in json_dict["superfeatures"]:
        assert list(superfeature.keys()) == superfeature_keys
        assert isinstance(superfeature["atom_numbers"], list)
        assert isinstance(superfeature["envpartners"], list)
        assert isinstance(superfeature["feature_type"], str)
        assert isinstance(superfeature["id"], str)
        assert isinstance(superfeature["occurrences"], list)

        for envpartner in superfeature["envpartners"]:
            assert list(envpartner.keys()) == envpartner_keys
            assert isinstance(envpartner["atom_numbers"], list)
            assert isinstance(envpartner["distances"], list)
            assert isinstance(envpartner["id"], str)
            assert isinstance(envpartner["name"], str)
            assert isinstance(envpartner["occurrences"], list)

    assert isinstance(json_dict["ligand_name"], str)
    assert isinstance(json_dict["ligand_smiles"], str)
