"""
Unit tests for dynophore.core.ligand.Ligand class.

Uses fixture tests.conftest.envpartner.
"""


from pathlib import Path

import pytest
from rdkit import Chem
from IPython.display import Image

from dynophores.core.ligand import Ligand

PATH_TEST_DATA = Path(__name__).parent / "dynophores/tests/data"


class TestsLigand:
    """
    Test Ligand class methods.
    """

    def test_from_dynophore(self, ligand, dynophore):

        assert list(ligand.__dict__) == [
            "name",
            "smiles",
            "_superfeatures_pdb_atom_ids",
            "_superfeatures_colors",
            "_pdb_block",
            "_pdb_atom_ids",
        ]
        assert isinstance(ligand, Ligand)
        assert isinstance(ligand.name, str)
        assert len(ligand.name) == 3
        assert isinstance(ligand.smiles, str)

        pdb_atom_ids = ligand._superfeatures_pdb_atom_ids
        assert isinstance(pdb_atom_ids, dict)
        assert list(pdb_atom_ids) == list(dynophore.superfeatures)
        assert isinstance(next(iter(pdb_atom_ids.values()))[0], int)

        colors = ligand._superfeatures_colors
        assert isinstance(colors, dict)
        assert list(colors) == list(dynophore.superfeatures)
        assert isinstance(next(iter(colors.values())), str)
        assert len(next(iter(colors.values()))) in [6, 7]  # With/without #-prefix

        assert isinstance(ligand._pdb_block, str)
        assert isinstance(ligand._pdb_atom_ids[0], int)

    @pytest.mark.parametrize("show_pdb_serial_numbers", [False, True])
    def test_view2d(self, ligand, show_pdb_serial_numbers):

        mol = ligand.view2d(show_pdb_serial_numbers)
        assert isinstance(mol, Chem.rdchem.Mol)

    @pytest.mark.parametrize("show_pdb_serial_numbers", [False, True])
    def test_view2d_superfeatures(self, ligand, show_pdb_serial_numbers):

        img = ligand.view2d_superfeatures(show_pdb_serial_numbers)
        assert isinstance(img, Image)
