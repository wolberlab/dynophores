"""
Unit tests for dynophore.viz.view2d.interactive.

Will only test if function signature is correct.
Does not test if called function executes without error (but those functions are tested elsewhere
anyways).
Issue reported here: https://github.com/jupyter-widgets/ipywidgets/issues/2949
"""

from pathlib import Path

import pytest

from dynophores.viz import view2d

PATH_TEST_DATA = Path(__name__).parent / "dynophores/tests/data"


@pytest.mark.parametrize("pdb_path", [PATH_TEST_DATA / "in/startframe.pdb"])
def test_show(dynophore, pdb_path):
    view2d.interactive.show(dynophore, pdb_path)
