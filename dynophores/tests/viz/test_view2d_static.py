"""
Unit tests for dynophore.viz.view2d.static.

Will only test if static 2D view raises errors.
"""

from pathlib import Path

import pytest
from rdkit import Chem
from IPython.display import Image

from dynophores.viz import view2d

PATH_TEST_DATA = Path(__name__).parent / "dynophores/tests/data"


@pytest.mark.parametrize(
    "show_superfeatures, show_pdb_serial_numbers",
    [
        (False, False),
        (False, True),
        (True, False),
        (True, True),
    ],
)
def test_show(dynophore, show_superfeatures, show_pdb_serial_numbers):
    view2d.static.show(dynophore, show_superfeatures, show_pdb_serial_numbers)

    @pytest.mark.parametrize("show_pdb_serial_numbers", [False, True])
    def test_view2d(self, ligand, show_pdb_serial_numbers):
        mol = ligand.view2d(show_pdb_serial_numbers)
        assert isinstance(mol, Chem.rdchem.Mol)

    @pytest.mark.parametrize("show_pdb_serial_numbers", [False, True])
    def test_view2d_superfeatures(self, ligand, show_pdb_serial_numbers):
        img = ligand.view2d_superfeatures(show_pdb_serial_numbers)
        assert isinstance(img, Image)
