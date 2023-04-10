"""
Unit tests for dynophore.core.ligand.Ligand class.

Uses fixture tests.conftest.envpartner.
"""


from pathlib import Path

from dynophores import parsers
from dynophores.core.ligand import Ligand

PATH_TEST_DATA = Path(__name__).parent / "dynophores/tests/data"


class TestsLigand:
    """
    Test Ligand class methods.
    """

    def test_init(self):
        dynophore_dict = parsers._json_pml_to_dict(
            PATH_TEST_DATA / "out/1KE7_dynophore.json",
            PATH_TEST_DATA / "out/1KE7_dynophore.pml",
        )
        ligand_dict = dynophore_dict["ligand"]
        ligand = Ligand(**ligand_dict)

        assert list(ligand.__dict__) == [
            "name",
            "smiles",
            "mdl_mol_buffer",
            "atom_serials",
        ]
        assert isinstance(ligand, Ligand)
        assert isinstance(ligand.name, str)
        assert len(ligand.name) == 3
        assert isinstance(ligand.smiles, str)
        assert isinstance(ligand.mdl_mol_buffer, str)
        assert isinstance(ligand.atom_serials[0], int)
