"""
Python API to create (i.e. generate and visualize) dynophore data.
"""

from pathlib import Path
import subprocess
import glob

from ..utils import _download_file
from ..definitions import DYNOPHORE_JAR_URL, DYNOPHORE_JAR_PATH
from .visualize import visualize


def create(pmz_or_pdb_path, dcd_path, out_path, name, feature_def_path, three_letter_code, chain):
    """
    Generate and visualize (Jupyter notebook) dynophore data.

    Parameters
    ----------
    pmz_or_pdb_path : str or pathlib.Path
        PMZ or PDB file.
    dcd_path : str or pathlib.Path
        DCD file.
    out_path : str or pathlib.Path
        Output folder.
    name : str
        Dynophore name.
    feature_def_path : str
        Optional: Feature definition file.
    three_letter_code : str
        Optional: Three letter code of ligand to be used for binding site definition.
        In combination with PDB file.
    chain : str
        Optional: Chain to be used for binding site definition.
        In combination with PDB file.
    """

    # Generate dynophore data
    generate(pmz_or_pdb_path, dcd_path, out_path, name, feature_def_path, three_letter_code, chain)

    # Visualize dynophore data
    if pmz_or_pdb_path.suffix == ".pdb":
        pdb_path = pmz_or_pdb_path
        dyno_paths = sorted(glob.glob(str(out_path / "dynophore_out_*")), reverse=True)
        dyno_path = Path(dyno_paths[0])
        visualize(dyno_path, pdb_path, dcd_path)
    else:
        print("To generate viz, provide PDB file.")


def generate(
    pmz_or_pdb_path,
    dcd_path,
    out_path,
    name,
    feature_def_path=None,
    three_letter_code=None,
    chain=None,
):
    """
    Generate dynophore data using the dynophore Java library.

    Parameters
    ----------
    pmz_or_pdb_path : str or pathlib.Path
        PMZ or PDB file.
    dcd_path : str or pathlib.Path
        DCD file.
    out_path : str or pathlib.Path
        Output folder.
    name : str
        Dynophore name.
    feature_def_path : str
        Optional: Feature definition file.
    three_letter_code : str
        Optional: Three letter code of ligand to be used for binding site definition.
        In combination with PDB file.
    chain : str
        Optional: Chain to be used for binding site definition.
        In combination with PDB file.
    """

    _download_file(DYNOPHORE_JAR_URL, DYNOPHORE_JAR_PATH)

    command = (
        f"java -jar {DYNOPHORE_JAR_PATH} "
        f"--pmz {pmz_or_pdb_path} "
        f"--dcd {dcd_path} "
        f"--out {out_path} "
        f"--name {name} "
        f"{'' if feature_def_path is None else f'--feature-def-file {feature_def_path}'} "
        f"{'' if three_letter_code is None else f'--three-letter-code {three_letter_code}'} "
        f"{'' if chain is None else f'--chain {chain}'} "
    )
    command = command.split()
    subprocess.run(command)
