"""
Python API to visualize dynophore data in a Jupyter notebook.
"""

from pathlib import Path

from .utils import _copy_notebook, _update_paths_in_notebook, _open_notebook

PATH_TEST_DATA = Path(__file__).parent / "tests" / "data"


def visualize(dyno_path, pdb_path, dcd_path=None):
    """
    Visualize dynophore data in a Jupyter notebook,

    Parameters
    ----------
    dyno_path : str or pathlib.Path
        Dynophore output folder.
    pdb_path : str or pathlib.Path
        PDB file.
    dcd_path : str or pathlib.Path
        Optional: DCD file.
    """

    new_notebook_path = dyno_path / "dynophore.ipynb"

    if not dyno_path.is_dir():
        raise RuntimeError(
            f"Input is no directory or directory does not exist: `{dyno_path.absolute()}`"
        )

    if not pdb_path.is_file():
        raise RuntimeError(f"Input is no file or file does not exist: `{pdb_path.absolute()}`")

    if dcd_path is not None:
        if not dcd_path.is_file():
            raise RuntimeError(f"Input is no file or file does not exist: `{dcd_path.absolute()}`")

    _copy_notebook(new_notebook_path)
    _update_paths_in_notebook(new_notebook_path, dyno_path, pdb_path, dcd_path)
    _open_notebook(new_notebook_path)
