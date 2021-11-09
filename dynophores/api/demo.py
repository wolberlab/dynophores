"""
"""

import json
from pathlib import Path
from shutil import copyfile
import subprocess


from .. import _version
from .utils import _copy_notebook, _update_paths_in_notebook, _open_notebook

PATH_TEST_DATA = Path(__file__).parent / ".." / "tests" / "data"


def demo(out_path):
    """
    Create and open demo visualization notebook based.

    Parameters
    ----------
    out_path : str or pathlib.Path
        Output folder.
    """

    new_notebook_path = Path(out_path) / "dynophore_demo.ipynb"
    _copy_notebook(new_notebook_path)
    _update_paths_in_notebook(
        new_notebook_path,
        (PATH_TEST_DATA / "out").absolute(),
        (PATH_TEST_DATA / "in/startframe.pdb").absolute(),
        (PATH_TEST_DATA / "in/trajectory.dcd").absolute(),
    )
    _open_notebook(new_notebook_path)
