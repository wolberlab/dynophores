"""
Unit tests for the dynophore CLI.
"""

from pathlib import Path
import subprocess
import warnings
import glob
import shutil

import pytest

PATH_TEST_DATA = Path(__name__).parent / "dynophores" / "tests" / "data"


def capture(command):
    """
    Run input command as subprocess and caputure the subprocess' exit code, stdout and stderr.
    Parameters
    ----------
    command : list of str
        Command to be run as subprocess.
    Returns
    -------
    out : str
        Standard output message.
    err : str
        Standard error message.
    exitcode : int
        Exit code.
    """

    proc = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    out, err = proc.communicate()
    return out, err, proc.returncode


@pytest.mark.parametrize(
    "command",
    [
        ["dynophore", "demo", PATH_TEST_DATA],
    ],
)
def test_cli_demo(command):
    """
    Test CLI.
    """

    _, err, exitcode = capture(command)
    assert exitcode == 0
    assert not err

    # Clean up output
    (PATH_TEST_DATA / "dynophore_demo.ipynb").unlink()


@pytest.mark.parametrize(
    "command",
    [
        [  # Create output from PMZ
            "dynophore",
            "create",
            "-p",
            PATH_TEST_DATA / "in/startframe.pmz",
            "-d",
            PATH_TEST_DATA / "in/trajectory.dcd",
            "-o",
            PATH_TEST_DATA / "out",
            "-n",
            "test",
        ],
        [  # Create output from PDB
            "dynophore",
            "create",
            "-p",
            PATH_TEST_DATA / "in/startframe.pdb",
            "-d",
            PATH_TEST_DATA / "in/trajectory.dcd",
            "-o",
            PATH_TEST_DATA / "out",
            "-n",
            "test",
        ],
        [  # Create output from PDB with user-defined chain and ligand
            "dynophore",
            "create",
            "-p",
            PATH_TEST_DATA / "in/startframe.pdb",
            "-d",
            PATH_TEST_DATA / "in/trajectory.dcd",
            "-o",
            PATH_TEST_DATA / "out",
            "-n",
            "test",
            "-c",
            "A",
            "-3",
            "LS3",
        ],
        # TODO feature definitions
    ],
)
def test_cli_create(command):
    """
    Test CLI.
    """

    _, err, exitcode = capture(command)

    assert exitcode == 0

    # assert not err
    # FIXME at the moment there is a warning in stderr
    # Uncomment lines above once this is fixed
    assert (
        err
        == "jar:file:/home/dominique/Documents/GitHub/dynophores/dynophores/generate/dynophore-20201007.jar!/com/inteligand/ilib/fonts/mini.flf\n"
    )

    # Clean up output
    dyno_paths = sorted(glob.glob(str(PATH_TEST_DATA / "out" / "dynophore_out_*")), reverse=True)
    dyno_path = Path(dyno_paths[0])
    shutil.rmtree(dyno_path)


@pytest.mark.parametrize(
    "command",
    [
        [  # Visualize dynophore with PDB _and_ DCD
            "dynophore",
            "visualize",
            "-i",
            PATH_TEST_DATA / "out",
            "-p",
            PATH_TEST_DATA / "in/startframe.pdb",
            "-d",
            PATH_TEST_DATA / "in/trajectory.dcd",
        ],
    ],
)
def test_cli_visualize(command):
    """
    Test CLI.
    """

    _, err, exitcode = capture(command)
    assert exitcode == 0
    assert not err

    # Clean up output
    (PATH_TEST_DATA / "out" / "dynophore.ipynb").unlink()


@pytest.mark.parametrize(
    "command, out_substring",
    [
        (
            [  # Input PMZ unknown
                "dynophore",
                "create",
                "-p",
                PATH_TEST_DATA / "in/xxx.pmz",
                "-d",
                PATH_TEST_DATA / "in/trajectory.dcd",
                "-o",
                PATH_TEST_DATA / "out",
                "-n",
                "test",
            ],
            "file not found",
        ),
        (
            [  # Input PDB unknown
                "dynophore",
                "create",
                "-p",
                PATH_TEST_DATA / "in/xxx.pdb",
                "-d",
                PATH_TEST_DATA / "in/trajectory.dcd",
                "-o",
                PATH_TEST_DATA / "out",
                "-n",
                "test",
            ],
            "file not found",
        ),
        (
            [  # Input DCD unknown
                "dynophore",
                "create",
                "-p",
                PATH_TEST_DATA / "in/startframe.pmz",
                "-d",
                PATH_TEST_DATA / "in/xxx.dcd",
                "-o",
                PATH_TEST_DATA / "out",
                "-n",
                "test",
            ],
            "file not found",
        ),
        (
            [  # Input chain unknown
                "dynophore",
                "create",
                "-p",
                PATH_TEST_DATA / "in/startframe.pdb",
                "-d",
                PATH_TEST_DATA / "in/trajectory.dcd",
                "-o",
                PATH_TEST_DATA / "out",
                "-n",
                "test",
                "-c",
                "X",
            ],
            "ERROR: No contained ligand matches given specification",
        ),
        (
            [  # Input ligand unknown
                "dynophore",
                "create",
                "-p",
                PATH_TEST_DATA / "in/startframe.pdb",
                "-d",
                PATH_TEST_DATA / "in/trajectory.dcd",
                "-o",
                PATH_TEST_DATA / "out",
                "-n",
                "test",
                "-3",
                "XXX",
            ],
            "ERROR: No contained ligand matches given specification",
        ),
        # TODO feature definitions
    ],
)
def test_cli_create_raises(command, out_substring):
    """
    Test CLI errors.
    """

    out, err, exitcode = capture(command)
    # Since the dynophore is generated running a JAR file,
    # we won't see errors from our Python program
    assert exitcode == 0

    # assert not err
    # FIXME at the moment there is a warning in stderr
    # Uncomment lines above once this is fixed
    assert (
        err
        == "jar:file:/home/dominique/Documents/GitHub/dynophores/dynophores/generate/dynophore-20201007.jar!/com/inteligand/ilib/fonts/mini.flf\n"
    )
    assert out_substring in out

    # Clean up output
    dyno_paths = sorted(glob.glob(str(PATH_TEST_DATA / "out" / "dynophore_out_*")), reverse=True)
    dyno_path = Path(dyno_paths[0])
    shutil.rmtree(dyno_path)


@pytest.mark.parametrize(
    "command",
    [
        [  # Dynophore folder unknown
            "dynophore",
            "visualize",
            "-i",
            PATH_TEST_DATA / "xxx",
            "-p",
            PATH_TEST_DATA / "in/startframe.pdb",
        ],
        [  # Input PDB unknown
            "dynophore",
            "visualize",
            "-i",
            PATH_TEST_DATA / "out",
            "-p",
            PATH_TEST_DATA / "in/xxx.pdb",
        ],
        [  # Input DCD unknown
            "dynophore",
            "visualize",
            "-i",
            PATH_TEST_DATA / "out",
            "-p",
            PATH_TEST_DATA / "in/startframe.pdb",
            "-d",
            PATH_TEST_DATA / "in/xxx.dcd",
        ],
    ],
)
def test_cli_visualize_raises(command):
    """
    Test CLI errors.
    """

    _, _, exitcode = capture(command)
    assert exitcode == 1


def test_jlab_import():
    """
    Add warning if JupyterLab cannot be imported.
    """

    try:
        import jupyterlab
    except ImportError:
        warnings.warn(
            "JupyterLab cannot be imported; install with `mamba install jupyterlab`.",
            ImportWarning,
        )
