"""
Unit tests for the dynophore CLI.
"""

from pathlib import Path
import subprocess

import pytest

from dynophores import cli

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

