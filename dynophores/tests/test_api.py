"""
Unit tests for the dynophore API.
"""

from pathlib import Path
import glob
import shutil

import pytest

from dynophores import api

PATH_TEST_DATA = Path(__name__).parent / "dynophores" / "tests" / "data"


@pytest.mark.parametrize(
    "pmz_path, dcd_path, out_path, name, feature_def_path, three_letter_code, chain",
    [
        (
            PATH_TEST_DATA / "in/startframe.pmz",
            PATH_TEST_DATA / "in/trajectory.dcd",
            PATH_TEST_DATA / "out",
            "test",
            None,
            None,
            None,
        ),
    ],
)
def test_create_from_pmz(
    pmz_path, dcd_path, out_path, name, feature_def_path, three_letter_code, chain
):
    """
    Test creating dynophore data from a PMZ file.
    """

    api.create.create(
        pmz_path, dcd_path, out_path, name, feature_def_path, three_letter_code, chain
    )

    # Check output files
    dyno_paths = sorted(glob.glob(str(out_path / "dynophore_out_*")), reverse=True)
    dyno_path = Path(dyno_paths[0])
    assert (dyno_path / f"{name}_dynophore.cgo").exists()
    assert (dyno_path / f"{name}_dynophore.json").exists()
    assert (dyno_path / f"{name}_dynophore.pml").exists()
    assert (dyno_path / f"{name}_input_filepaths.txt").exists()
    assert not (dyno_path / f"dynophore.ipynb").exists()
    shutil.rmtree(dyno_path)


@pytest.mark.parametrize(
    "pdb_path, dcd_path, out_path, name, feature_def_path, three_letter_code, chain",
    [
        (
            PATH_TEST_DATA / "in/startframe.pdb",
            PATH_TEST_DATA / "in/trajectory.dcd",
            PATH_TEST_DATA / "out",
            "test",
            None,
            None,
            None,
        ),
        (
            PATH_TEST_DATA / "in/startframe.pmz",
            PATH_TEST_DATA / "in/trajectory.dcd",
            PATH_TEST_DATA / "out",
            "test",
            None,
            "LS3",
            "A",
        ),
    ],
)
def test_create_from_pdb_xxx(
    pdb_path, dcd_path, out_path, name, feature_def_path, three_letter_code, chain
):
    """
    Test creating dynophore data from a PDB file.
    """

    api.create.create(
        pdb_path, dcd_path, out_path, name, feature_def_path, three_letter_code, chain
    )

    # Check output files
    dyno_paths = sorted(glob.glob(str(out_path / "dynophore_out_*")), reverse=True)
    dyno_path = Path(dyno_paths[0])
    assert (dyno_path / f"{name}_dynophore.cgo").exists()
    assert (dyno_path / f"{name}_dynophore.json").exists()
    assert (dyno_path / f"{name}_dynophore.pml").exists()
    assert (dyno_path / f"{name}_input_filepaths.txt").exists()
    assert (dyno_path / f"dynophore.ipynb").exists()
    shutil.rmtree(dyno_path)


@pytest.mark.parametrize(
    "pdb_path, dcd_path, out_path, name, feature_def_path, three_letter_code, chain",
    [
        (
            PATH_TEST_DATA / "in/startframe.pdb",
            PATH_TEST_DATA / "in/trajectory.dcd",
            PATH_TEST_DATA / "out",
            "test",
            None,
            "XXX",
            None,
        ),
        (
            PATH_TEST_DATA / "in/startframe.pmz",
            PATH_TEST_DATA / "in/trajectory.dcd",
            PATH_TEST_DATA / "out",
            "test",
            None,
            None,
            "X",
        ),
    ],
)
def test_create_from_pdb_invalid(
    pdb_path, dcd_path, out_path, name, feature_def_path, three_letter_code, chain
):
    """
    Test creating dynophore data from a PDB file with invalid user-input 3-letter code or chain.
    """

    api.create.create(
        pdb_path, dcd_path, out_path, name, feature_def_path, three_letter_code, chain
    )

    # Check output files
    dyno_paths = sorted(glob.glob(str(out_path / "dynophore_out_*")), reverse=True)
    dyno_path = Path(dyno_paths[0])
    assert not (dyno_path / f"{name}_dynophore.cgo").exists()
    assert not (dyno_path / f"{name}_dynophore.json").exists()
    assert not (dyno_path / f"{name}_dynophore.pml").exists()
    assert not (dyno_path / f"{name}_input_filepaths.txt").exists()
    assert not (dyno_path / f"dynophore.ipynb").exists()
    shutil.rmtree(dyno_path)


@pytest.mark.parametrize(
    "dyno_path, pdb_path, dcd_path",
    [
        (
            PATH_TEST_DATA / "out",
            PATH_TEST_DATA / "in/startframe.pdb",
            PATH_TEST_DATA / "in/trajectory.dcd",
        ),
        (
            PATH_TEST_DATA / "out",
            PATH_TEST_DATA / "in/startframe.pdb",
            None,
        ),
    ],
)
def test_visualize(dyno_path, pdb_path, dcd_path):
    """
    Test visualizing dynophore data.
    """

    api.visualize.visualize(dyno_path, pdb_path, dcd_path)
    assert (dyno_path / f"dynophore.ipynb").exists()
    (dyno_path / "dynophore.ipynb").unlink()


@pytest.mark.parametrize(
    "workspace_path",
    [(PATH_TEST_DATA)],
)
def test_demo(workspace_path):
    """
    Test demo notebook generation.
    """

    api.demo.demo(workspace_path)
    assert (workspace_path / f"dynophore_demo.ipynb").exists()
    (workspace_path / "dynophore_demo.ipynb").unlink()


@pytest.mark.parametrize(
    "new_notebook_path",
    [PATH_TEST_DATA / "copied_notebook.ipynb"],
)
def test_copy_notebook(new_notebook_path):
    """
    Test if copied notebook file exists.
    """

    api.utils._copy_notebook(new_notebook_path)
    assert new_notebook_path.is_file()
    # Remove copy again (not needed)
    new_notebook_path.unlink()


@pytest.mark.parametrize(
    "new_notebook_path",
    ["xxx"],
)
def test_copy_notebook_raises(new_notebook_path):
    """
    Test if error raised if input path does not have suffix ".ipynb".
    """

    with pytest.raises(RuntimeError):
        api.utils._copy_notebook(new_notebook_path)


@pytest.mark.parametrize(
    "notebook_path, dyno_path, pdb_path, dcd_path",
    [
        (
            str(PATH_TEST_DATA / "test.ipynb"),  # Not a file
            str(PATH_TEST_DATA / "out"),
            str(PATH_TEST_DATA / "in/startframe.pdb"),
            str(PATH_TEST_DATA / "in/trajectory.dcd"),
        ),
        (
            str(PATH_TEST_DATA / "test.ipynb"),  # TODO
            "is_not_dir",  # Not a directory
            str(PATH_TEST_DATA / "in/startframe.pdb"),
            str(PATH_TEST_DATA / "in/trajectory.dcd"),
        ),
        (
            str(PATH_TEST_DATA / "test.ipynb"),  # TODO
            str(PATH_TEST_DATA / "out"),
            "doesnt_exist.pdb",  # Does not exist
            str(PATH_TEST_DATA / "in/trajectory.dcd"),
        ),
        (
            str(PATH_TEST_DATA / "test.ipynb"),
            str(PATH_TEST_DATA / "out"),
            str(PATH_TEST_DATA / "in/startframe.pdb"),
            "doesnt_exist.dcd",  # Does not exist
        ),
    ],
)
def test_update_paths_in_notebook_raises(notebook_path, dyno_path, pdb_path, dcd_path):
    """
    Test if error is raised if input paths do not exist.
    """

    with pytest.raises(RuntimeError):
        api.utils._update_paths_in_notebook(notebook_path, dyno_path, pdb_path, dcd_path)
