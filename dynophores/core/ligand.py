"""
Handles the Ligand class, which contains information on a macromolecule-bound ligand.
"""

from rdkit import Chem
from rdkit.Chem import AllChem


class Ligand:
    """
    Class to store information on a macromolecule-bound ligand.

    Attributes
    ----------
    name : str
        3-letter name of structure-bound ligand.
    smiles : str
        SMILES of structure-bound ligand.
    mdl_mol_buffer : str
        Ligand PDB block
    atom_serials : list of int
        Ligand PDB atom IDs.
    """

    def __init__(
        self,
        name,
        smiles,
        mdl_mol_buffer,
        atom_serials,
        **kwargs,
    ):

        self.name = name
        self.smiles = smiles
        self.mdl_mol_buffer = mdl_mol_buffer
        self.atom_serials = atom_serials

    def rdkit_molecule(self, show_pdb_serial_numbers=False):
        """
        2D view of ligand. Uses rdkit.

        Parameters
        ----------
        show_pdb_serial_numbers : bool
            Show PDB serial numbers of ligand atoms (default: False).

        Returns
        -------
        mol : rdkit.Chem.rdchem.Mol
            Molecule.
        """

        # Read ligand from MDL mol buffer; compute 2D coordinates
        mol = Chem.rdmolfiles.MolFromMolBlock(self.mdl_mol_buffer, removeHs=False)
        AllChem.Compute2DCoords(mol)

        # Set PDB atom serial numbers to rdkit molecule
        if len(self.atom_serials) == mol.GetNumAtoms():
            for atom, pdb_atom_serial in zip(mol.GetAtoms(), self.atom_serials):
                atom.SetIntProp("pdb_atom_serial", int(pdb_atom_serial))
        else:
            raise ValueError(
                f"Number of atoms serials ({len(self.atom_serials)}) and "
                f"number of atoms in RDKit molecule ({mol.GetNumAtoms()}) are not the same."
            )

        # Remove hydrogens
        mol = Chem.RemoveHs(mol)

        if show_pdb_serial_numbers:
            for atom in mol.GetAtoms():
                atom.SetAtomMapNum(atom.GetIntProp("pdb_atom_serial"))

        return mol
