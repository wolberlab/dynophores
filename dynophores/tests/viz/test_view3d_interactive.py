"""
Unit tests for dynophore.viz.view3d.interactive.

Will only test if function signature is correct.
Does not test if called function executes without error (but those functions are tested elsewhere
anyways).
Issue reported here: https://github.com/jupyter-widgets/ipywidgets/issues/2949

"""

from pathlib import Path

import pytest

from dynophores.viz import view3d

PATH_TEST_DATA = Path(__name__).parent / "dynophores" / "tests" / "data"


@pytest.mark.parametrize(
    "pdb_path, dcd_path, visualization_type, color_cloud_by_frame",
    [
        (PATH_TEST_DATA / "in/startframe.pdb", None, "spheres", False),
        (
            PATH_TEST_DATA / "in/startframe.pdb",
            PATH_TEST_DATA / "in/trajectory.dcd",
            "points",
            False,
        ),
        (PATH_TEST_DATA / "in/startframe.pdb", None, "spheres", True),
    ],
)
def test_show(dynophore, pdb_path, dcd_path, visualization_type, color_cloud_by_frame):
    view3d.show(dynophore, pdb_path, dcd_path, visualization_type, color_cloud_by_frame)


@pytest.mark.parametrize(
    "pdb_path, dcd_path, visualization_type",
    [
        (PATH_TEST_DATA / "in/startframe.pdb", None, "xxx"),
    ],
)
def test_show_raises(dynophore, pdb_path, dcd_path, visualization_type):

    with pytest.raises(ValueError):
        view3d.show(dynophore, pdb_path, dcd_path, visualization_type)
