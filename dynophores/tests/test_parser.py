"""
Unit tests for dynophore.viz.view3d.parser.
"""

from pathlib import Path

import pytest
import numpy as np

from dynophores import parsers

PATH_TEST_DATA = Path(__name__).parent / "dynophores" / "tests" / "data"


@pytest.mark.parametrize(
    "filepath, superfeature_names, cloud_keys",
    [
        (
            PATH_TEST_DATA / "out/1KE7_dynophore.pml",
            [
                "AR[4605,4607,4603,4606,4604]",
                "AR[4622,4615,4623,4613,4614,4621]",
                "HBA[4596]",
                "HBA[4606]",
                "HBA[4618]",
                "HBA[4619]",
                "HBD[4598]",
                "HBD[4612]",
                "H[4599,4602,4601,4608,4609,4600]",
                "H[4615,4623,4622,4613,4621,4614]",
            ],
            ["center", "coordinates", "id"],
        )
    ],
)
def test_pml_to_dict(filepath, superfeature_names, cloud_keys):

    pml_dict = parsers._pml_to_dict(filepath)
    assert sorted(pml_dict) == superfeature_names
    for superfeature_name, data in pml_dict.items():
        assert sorted(data.keys()) == cloud_keys
        assert isinstance(superfeature_name, str)
        assert data["id"] == superfeature_name
        assert isinstance(data["center"], np.ndarray)
        assert isinstance(data["coordinates"], np.ndarray)


@pytest.mark.parametrize(
    "filepath, dynophore_keys, superfeature_names, superfeature_keys, envpartner_keys",
    [
        (
            PATH_TEST_DATA / "out/1KE7_dynophore.json",
            ["id", "superfeatures"],
            [
                "AR[4605,4607,4603,4606,4604]",
                "AR[4622,4615,4623,4613,4614,4621]",
                "HBA[4596]",
                "HBA[4606]",
                "HBA[4618]",
                "HBA[4619]",
                "HBD[4598]",
                "HBD[4612]",
                "H[4599,4602,4601,4608,4609,4600]",
                "H[4615,4623,4622,4613,4621,4614]",
            ],
            [
                "atom_numbers",
                "envpartners",
                "feature_type",
                "id",
                "occurrences",
            ],
            ["atom_numbers", "distances", "id", "name", "occurrences"],
        )
    ],
)
def test_json_to_dict(
    filepath, dynophore_keys, superfeature_names, superfeature_keys, envpartner_keys
):

    pml_dict = parsers._json_to_dict(filepath)
    assert sorted(pml_dict.keys()) == dynophore_keys
    assert isinstance(pml_dict["superfeatures"], list)

    for superfeature in pml_dict["superfeatures"]:
        assert sorted(superfeature.keys()) == superfeature_keys
        assert isinstance(superfeature["atom_numbers"], list)
        assert isinstance(superfeature["envpartners"], list)
        assert isinstance(superfeature["feature_type"], str)
        assert isinstance(superfeature["id"], str)
        assert isinstance(superfeature["occurrences"], list)

        for envpartner in superfeature["envpartners"]:
            assert sorted(envpartner.keys()) == envpartner_keys
            assert isinstance(envpartner["atom_numbers"], list)
            assert isinstance(envpartner["distances"], list)
            assert isinstance(envpartner["id"], str)
            assert isinstance(envpartner["name"], str)
            assert isinstance(envpartner["occurrences"], list)
