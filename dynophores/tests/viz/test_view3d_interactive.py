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
    "pdb_path, dcd_path, visualization_type, color_cloud_by_frame, frame_range",
    [
        (PATH_TEST_DATA / "in/startframe.pdb", None, "spheres", False, None),
        (
            PATH_TEST_DATA / "in/startframe.pdb",
            PATH_TEST_DATA / "in/trajectory.dcd",
            "points",
            False,
            None,
        ),
        (PATH_TEST_DATA / "in/startframe.pdb", None, "spheres", True, None),
        (PATH_TEST_DATA / "in/startframe.pdb", None, "spheres", False, [100, 200]),
    ],
)
def test_show(
    dynophore, pdb_path, dcd_path, visualization_type, color_cloud_by_frame, frame_range
):
    view3d.show(
        dynophore,
        pdb_path,
        dcd_path=dcd_path,
        visualization_type=visualization_type,
        color_cloud_by_frame=color_cloud_by_frame,
        frame_range=frame_range,
    )


@pytest.mark.parametrize(
    "pdb_path, dcd_path, visualization_type, color_cloud_by_frame, frame_range",
    [
        (PATH_TEST_DATA / "in/startframe.pdb", None, "xxx", None, None),  # Unknown viz type
        (
            PATH_TEST_DATA / "in/startframe.pdb",
            None,
            "spheres",
            None,
            [0],
        ),  # Range needs two numbers
        (
            PATH_TEST_DATA / "in/startframe.pdb",
            None,
            "spheres",
            None,
            [100, 0],
        ),  # Range needs increasing numbers
    ],
)
def test_show_raises(
    dynophore, pdb_path, dcd_path, visualization_type, color_cloud_by_frame, frame_range
):
    with pytest.raises(ValueError):
        view3d.show(
            dynophore,
            pdb_path,
            dcd_path=dcd_path,
            visualization_type=visualization_type,
            color_cloud_by_frame=color_cloud_by_frame,
            frame_range=frame_range,
        )
