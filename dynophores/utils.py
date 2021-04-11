"""
Utility functions.
"""

import os
import shutil
import tempfile
import contextlib
from pathlib import Path

import MDAnalysis as mda


@contextlib.contextmanager
def enter_temp_directory(remove=True):
    """Create and enter a temporary directory; used as context manager."""
    temp_dir = tempfile.mkdtemp()
    cwd = os.getcwd()
    os.chdir(temp_dir)
    yield cwd, temp_dir
    os.chdir(cwd)
    if remove:
        shutil.rmtree(temp_dir)


def hex_to_rgb(hex_string, scale=True):
    """
    Convert hex code to RGB code.

    Parameters
    ----------
    hex_string : str
        6-letter hex string.

    Return
    ------
    list of int or float
        RGB code, e.g. (255, 255, 255).
    """

    if not isinstance(hex_string, str):
        raise ValueError(f"Input must be string but is {type(hex_string)}.")
    hex_string = hex_string.lstrip("#")
    if len(hex_string) != 6:
        raise ValueError(f"Input hex string must be of length 6 but is '{hex_string}'.")
    r_hex = hex_string[0:2]
    g_hex = hex_string[2:4]
    b_hex = hex_string[4:6]

    rgb = [int(r_hex, 16), int(g_hex, 16), int(b_hex, 16)]

    if scale:
        rgb = [round(value / float(255), 3) for value in rgb]

    return rgb


def pdb_ligand_data_for_rdkit(pdb_path, ligand_name):
    """
    Extract from structure PDB file ligand by ligand name.

    Parameters
    ----------
    pdb_path : str or pathlib.Path
        Path to PDB file
    ligand_name : str
        3-letter name of structure-bound ligand.

    Returns
    -------
    pdb_block : str
        Ligand PDB block
    pdb_atom_ids : list of int
        Ligand PDB atom IDs.
    """

    pdb_path = Path(pdb_path)

    # 1. Extract ligand from PDB using MDAnalysis
    universe = mda.Universe(pdb_path)
    ligand_mda = universe.select_atoms(f"resname {ligand_name}")
    if len(ligand_mda.atoms) == 0:
        raise ValueError("Input ligand name is not part of input PDB file.")

    # 2. Save ligand as PDB temporily to load it into `rdkit`
    ligand_pdb_atom_ids = [int(atom.id) for atom in ligand_mda.atoms]
    with enter_temp_directory():
        ligand_mda.write("ligand.pdb")
        with open("ligand.pdb", "r") as f:
            ligand_pdb_block = f.read()

    return ligand_pdb_block, ligand_pdb_atom_ids
