"""
Unit tests for Command Line Interface.
"""

from pathlib import Path
import subprocess

import pytest

from dynophores import cli

PATH_TEST_DATA = Path(__name__).parent / "dynophores" / "tests" / "data"


@pytest.mark.parametrize(
    "function, dyno, pdb, dcd, workspace",
    [
        (
            "create",
            PATH_TEST_DATA / "1KE7-1/DynophoreApp",
            PATH_TEST_DATA / "1KE7-1/startframe.pdb",
            PATH_TEST_DATA / "1KE7-1/trajectory.dcd",
            PATH_TEST_DATA,
        ),
    ],
)
def test_create_subprocess(function, dyno, pdb, dcd, workspace):
    """
    TODO how to catch errors?
    """

    try:
        subprocess.run(
            [
                "dynoviz",
                function,
                "--dyno",
                dyno,
                "--pdb",
                pdb,
                "--dcd",
                dcd,
                "--workspace",
                workspace,
            ],
            check=True,
            timeout=5,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )
    except subprocess.TimeoutExpired:
        pass


@pytest.mark.parametrize(
    "function, dyno, pdb, dcd, workspace, error",
    [
        ("create", "xxx", "xxx", "xxx", "xxx", subprocess.CalledProcessError),
        ("cook", "xxx", "xxx", "xxx", "xxx", subprocess.CalledProcessError),
    ],
)
def test_create_subprocess_raises(function, dyno, pdb, dcd, workspace, error):

    with pytest.raises(error):
        process = subprocess.run(
            [
                "dynoviz",
                function,
                "--dyno",
                dyno,
                "--pdb",
                pdb,
                "--dcd",
                dcd,
                "--workspace",
                workspace,
            ],
            check=True,
        )
        process.terminate()


@pytest.mark.parametrize(
    "function, notebook",
    [
        (
            "open",
            PATH_TEST_DATA / "dynophore.ipynb",
        )
    ],
)
def test_open_subprocess(function, notebook):

    process = subprocess.Popen(
        [
            "dynoviz",
            function,
            "notebook",
            notebook,
        ],
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )
    process.terminate()


@pytest.mark.parametrize(
    "workspace, dyno, pdb, dcd",
    [
        (
            "xxx",
            PATH_TEST_DATA / "1KE7-1/DynophoreApp",
            PATH_TEST_DATA / "1KE7-1/startframe.pdb",
            PATH_TEST_DATA / "1KE7-1/trajectory.dcd",
        ),
        (
            "xxx",
            "xxx",
            PATH_TEST_DATA / "1KE7-1/startframe.pdb",
            PATH_TEST_DATA / "1KE7-1/trajectory.dcd",
        ),
        (
            "xxx",
            "xxx",
            "xxx",
            PATH_TEST_DATA / "1KE7-1/trajectory.dcd",
        ),
        ("xxx", "xxx", "xxx", "xxx"),
    ],
)
def test_copy_notebook_raises(workspace, dyno, pdb, dcd):

    with pytest.raises(RuntimeError):
        cli._copy_notebook(workspace, dyno, pdb, dcd)


@pytest.mark.parametrize("notebook", [("xxx")])
def test_open_notebook_raises(notebook):

    with pytest.raises(RuntimeError):
        cli._open_notebook(notebook)
