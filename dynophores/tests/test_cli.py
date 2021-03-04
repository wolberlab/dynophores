"""
Unit tests for Command Line Interface.
"""

from pathlib import Path
import subprocess

import pytest

from dynophores import cli

PATH_TEST_DATA = Path(__name__).parent / "dynophores" / "tests" / "data"


@pytest.mark.parametrize("function", ["cook"])
def test_subprocess_raises(function):
    """
    Test incorrect subprocess name.
    """

    with pytest.raises(subprocess.CalledProcessError):
        process = subprocess.run(["dynoviz", function], check=True)
        process.terminate()


@pytest.mark.parametrize("function", ["demo"])
def test_subprocess_demo(function):
    """
    TODO how to catch errors?
    """

    process = subprocess.Popen(
        ["dynoviz", function], text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE
    )
    process.terminate()


@pytest.mark.parametrize(
    "function, args_dyno, args_pdb, args_dcd, args_workspace",
    [
        (
            "create",
            str(PATH_TEST_DATA / "out"),
            str(PATH_TEST_DATA / "in/startframe.pdb"),
            str(PATH_TEST_DATA / "in/trajectory.dcd"),
            str(PATH_TEST_DATA),
        ),
    ],
)
def test_subprocess_create(function, args_dyno, args_pdb, args_dcd, args_workspace):
    """
    TODO how to catch errors?
    """

    process = subprocess.Popen(
        [
            "dynoviz",
            function,
            "--dyno",
            args_dyno,
            "--pdb",
            args_pdb,
            "--dcd",
            args_dcd,
            "--workspace",
            args_workspace,
        ],
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )
    process.terminate()


@pytest.mark.parametrize(
    "function, args_dyno, args_pdb, args_dcd, args_workspace",
    [
        ("create", "xxx", "xxx", "xxx", "xxx"),
        (
            "create",
            str(PATH_TEST_DATA / "out"),
            "xxx",
            "xxx",
            "xxx",
        ),
        (
            "create",
            str(PATH_TEST_DATA / "out"),
            str(PATH_TEST_DATA / "in/startframe.pdb"),
            "xxx",
            "xxx",
        ),
        (
            "create",
            str(PATH_TEST_DATA / "out"),
            str(PATH_TEST_DATA / "in/startframe.pdb"),
            str(PATH_TEST_DATA / "in/trajectory.dcd"),
            "xxx",
        ),
    ],
)
def test_subprocess_create_raises(function, args_dyno, args_pdb, args_dcd, args_workspace):
    """
    Test invalid CLI args for create subprocess.
    """

    with pytest.raises(subprocess.CalledProcessError):
        process = subprocess.run(
            [
                "dynoviz",
                function,
                "--dyno",
                args_dyno,
                "--pdb",
                args_pdb,
                "--dcd",
                args_dcd,
                "--workspace",
                args_workspace,
            ],
            check=True,
        )
        process.terminate()


@pytest.mark.parametrize(
    "function, args_notebook",
    [
        (
            "open",
            str(PATH_TEST_DATA / "dynophore.ipynb"),
        )
    ],
)
def test_subprocess_open(function, args_notebook):
    """
    TODO how to catch errors?
    """

    process = subprocess.Popen(
        [
            "dynoviz",
            function,
            "notebook",
            args_notebook,
        ],
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )
    process.terminate()


@pytest.mark.parametrize(
    "function, args_notebook",
    [
        (
            "open",
            str(PATH_TEST_DATA / "xxx.ipynb"),
        )
    ],
)
def test_subprocess_open_raises(function, args_notebook):
    """
    Test invalid CLI args for open subprocess.
    """

    with pytest.raises(subprocess.CalledProcessError):
        process = subprocess.run(
            ["dynoviz", function, args_notebook],
            check=True,
        )
        process.terminate()


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
            "xxx",
            str(PATH_TEST_DATA / "out"),
            str(PATH_TEST_DATA / "in/startframe.pdb"),
            str(PATH_TEST_DATA / "in/trajectory.dcd"),
        ),
        (
            "xxx",
            "xxx",
            str(PATH_TEST_DATA / "in/startframe.pdb"),
            str(PATH_TEST_DATA / "in/trajectory.dcd"),
        ),
        (
            "xxx",
            "xxx",
            "xxx",
            str(PATH_TEST_DATA / "in/trajectory.dcd"),
        ),
        ("xxx", "xxx", "xxx", "xxx"),
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
