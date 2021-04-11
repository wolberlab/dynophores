"""
Unit tests for dynophore.viz.view3d.interactive.

Will only test if viewer raises errors.
"""

from pathlib import Path

import pytest

from dynophores.viz import view3d

PATH_TEST_DATA = Path(__name__).parent / "dynophores" / "tests" / "data"


@pytest.mark.parametrize(
    "pdb_path, dcd_path",
    [
        (
            PATH_TEST_DATA / "in/startframe.pdb",
            None,
        ),
        (
            PATH_TEST_DATA / "in/startframe.pdb",
            PATH_TEST_DATA / "in/trajectory.dcd",
        ),
    ],
)
def test_show(dynophore, pdb_path, dcd_path):
    view3d.show(dynophore, pdb_path, dcd_path)
