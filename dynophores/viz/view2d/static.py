"""
Contains RDKit ligand 2D visualizations (static).
"""

from dynophores.core.ligand import Ligand


def show(dynophore, pdb_path, show_superfeatures=False, show_pdb_serial_numbers=False):
    """
    Show static ligand 2D view.

    Parameters
    ----------
    pdb_path : str or pathlib.Path
        Path to PDB file
    dynophore : dynophores.Dynophore
        Dynophore object.

    Returns
    -------
    rdkit.Chem.rdchem.Mol or IPython.core.display.Image
        Molecule (with/without superfeatures).
    """

    ligand = Ligand.from_pdb(pdb_path, dynophore)
    if show_superfeatures:
        return ligand.view2d_superfeatures(show_pdb_serial_numbers)
    else:
        return ligand.view2d(show_pdb_serial_numbers)
