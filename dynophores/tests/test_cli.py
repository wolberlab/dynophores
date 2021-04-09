"""
Unit tests for Command Line Interface.
"""

from pathlib import Path
import subprocess

import pytest

from dynophores import cli

PATH_TEST_DATA = Path(__name__).parent / "dynophores" / "tests" / "data"


@pytest.mark.parametrize(
    "args",
    [
        "dynoviz demo",
        f"dynoviz create --dyno {PATH_TEST_DATA / 'out'} --workspace {PATH_TEST_DATA} "
        f"--pdb {PATH_TEST_DATA / 'in/startframe.pdb'} "
        f"--dcd {PATH_TEST_DATA / 'in/trajectory.dcd'}",
        f"dynoviz create --dyno {PATH_TEST_DATA / 'out'} --workspace {PATH_TEST_DATA} "
        f"--pdb {PATH_TEST_DATA / 'in/startframe.pdb'}",
        f"dynoviz open {PATH_TEST_DATA / 'dynophore.ipynb'}",
    ],
)
def test_subprocess(args):
    """
    TODO how to catch errors?
    """

    args = args.split()
    process = subprocess.Popen(
        args,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )
    process.terminate()


@pytest.mark.parametrize(
    "args",
    [
        "dynoviz cook",
        f"dynoviz create --dyno {PATH_TEST_DATA / 'out'} --workspace xxx --pdb xxx --dcd xxx",
        f"dynoviz create --dyno {PATH_TEST_DATA / 'out'} --workspace {PATH_TEST_DATA} "
        f"--pdb xxx --dcd xxx",
        f"dynoviz create --dyno {PATH_TEST_DATA / 'out'} --workspace {PATH_TEST_DATA} "
        f"--pdb {PATH_TEST_DATA / 'in/startframe.pdb'} --dcd xxx",
        f"dynoviz open {PATH_TEST_DATA / 'xxx.ipynb'}",
    ],
)
def test_subprocess_raises(args):
    """
    Test invalid CLI args for create subprocess.
    """

    args = args.split()
    with pytest.raises(subprocess.CalledProcessError):
        subprocess.run(args, check=True)


@pytest.mark.parametrize(
    "new_notebook_path",
    [PATH_TEST_DATA / "copied_notebook.ipynb"],
)
def test_copy_notebook(new_notebook_path):
    """
    Test if copied notebook file exists.
    """

    cli._copy_notebook(new_notebook_path)
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
        cli._copy_notebook(new_notebook_path)


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
        cli._update_paths_in_notebook(notebook_path, dyno_path, pdb_path, dcd_path)


@pytest.mark.parametrize(
    "notebook",
    [
        PATH_TEST_DATA / "xxx.ipynb",  # Path does not exist
        PATH_TEST_DATA,  # Path is not a file
        PATH_TEST_DATA / "dynophore.iiipynb",  # Incorrect file suffix
    ],
)
def test_open_notebook_raises(notebook):

    with pytest.raises(RuntimeError):
        cli._open_notebook(notebook)
