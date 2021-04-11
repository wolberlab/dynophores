"""
Unit tests for dynophore.viz.view2d.static.

Will only test if static 2D view raises errors.
"""

from pathlib import Path

import pytest

from dynophores.viz import view2d

PATH_TEST_DATA = Path(__name__).parent / "dynophores/tests/data"


@pytest.mark.parametrize(
    "pdb_path, show_superfeatures, show_pdb_serial_numbers",
    [
        (PATH_TEST_DATA / "in/startframe.pdb", False, False),
        (PATH_TEST_DATA / "in/startframe.pdb", False, True),
        (PATH_TEST_DATA / "in/startframe.pdb", True, False),
        (PATH_TEST_DATA / "in/startframe.pdb", True, True),
    ],
)
def test_show(dynophore, pdb_path, show_superfeatures, show_pdb_serial_numbers):
    view2d.static.show(dynophore, pdb_path, show_superfeatures, show_pdb_serial_numbers)
