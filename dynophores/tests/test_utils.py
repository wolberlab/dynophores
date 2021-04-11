"""
Unit tests for dynophores.utils.
"""

from pathlib import Path

import pytest

from dynophores import utils

PATH_TEST_DATA = Path(__name__).parent / "dynophores/tests/data"


@pytest.mark.parametrize(
    "hex_string, scale, rgb",
    [
        ("ffc20e", True, [1.0, 0.761, 0.055]),
        ("ffc20e", False, [255, 194, 14]),
        ("#ffc20e", False, [255, 194, 14]),
    ],
)
def test_hex_to_rgb(hex_string, scale, rgb):
    """
    Test hex to RGB conversion.
    """

    assert utils.hex_to_rgb(hex_string, scale) == rgb


@pytest.mark.parametrize(
    "hex_string",
    ["a", 0],
)
def test_hex_to_rgb_raises(hex_string):
    """
    Test hex to RGB conversion.
    """

    with pytest.raises(ValueError):
        assert utils.hex_to_rgb(hex_string)


@pytest.mark.parametrize(
    "pdb_path, ligand_name",
    [(PATH_TEST_DATA / "in/startframe.pdb", "XXX")],
)
def test_pdb_ligand_data_for_rdkit_raises(pdb_path, ligand_name):

    with pytest.raises(ValueError):
        utils.pdb_ligand_data_for_rdkit(pdb_path, ligand_name)
