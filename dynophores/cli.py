"""
Command Line Interface for the project.
"""

import argparse
from pathlib import Path
from shutil import copyfile
import subprocess

from . import _version


def main():
    """
    Main CLI function with the following signatures:

    dynoviz create
      --dyno path/to/dyno/folder
      --pdb path/to/pdb/file
      --dcd path/to/dcd/file
      --workspace path/to/workspace/folder

    dynoviz open
      --notebook path/to/existing/notebook
    """

    _greet()

    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers()

    create_subparser = subparsers.add_parser("create")
    open_subparser = subparsers.add_parser("open")

    # Arguments and function to be called for sub-command encode
    create_subparser.add_argument(
        "-i",
        "--dyno",
        type=str,
        help="Path to DynophoreApp output folder",
        required=True,
    )
    create_subparser.add_argument(
        "-p",
        "--pdb",
        type=str,
        help="Path to pdb (topology) file from trajectory, e.g. first frame",
        required=True,
    )
    create_subparser.add_argument(
        "-d",
        "--dcd",
        type=str,
        help="Path to dcd (trajectory) file",
        required=True,
    )
    create_subparser.add_argument(
        "-w",
        "--workspace",
        type=str,
        help="Path to workspace folder",
        required=True,
    )
    create_subparser.set_defaults(func=_create_viz)

    # Arguments and function to be called for sub-command compare
    open_subparser.add_argument(
        "notebook",
        help="Path to dynophore notebook file",
    )
    open_subparser.set_defaults(func=_open_viz)

    args = parser.parse_args()
    args.func(args)


def _greet():
    """
    Print CLI greeting.
    """

    logo_path = Path(_version.__file__).parent / "../docs/_static/stegosaurus.txt"
    with open(logo_path, "r", encoding="ascii") as f:
        print(f.read())
    title_str = f"Dynophores CLI {_version.get_versions()['version']}"
    print(f"\n{title_str:^64}")


def _create_viz(args):
    """
    Create visualization notebook based on command line arguments.
    """

    _copy_notebook(args.workspace, args.dyno, args.pdb, args.dcd)

    notebook_path = Path(args.workspace) / "dynophore.ipynb"
    _open_notebook(notebook_path)


def _open_viz(args):
    """
    Open visualization notebook based on command line arguments.
    """

    _open_notebook(args.notebook)


def _copy_notebook(workspace_path, dyno_path, pdb_path, dcd_path):
    """
    Copy template dynophore notebook to user-defined workspace and update filepaths set in notebook
    to user-defined filepaths.
    """

    workspace_path = Path(workspace_path)
    dyno_path = Path(dyno_path)
    pdb_path = Path(pdb_path)
    dcd_path = Path(dcd_path)

    if not workspace_path.exists():
        raise RuntimeError(f"Input file does not exist: `{workspace_path.absolute()}`")
    if not dyno_path.exists():
        raise RuntimeError(f"Input file does not exist: `{dyno_path.absolute()}`")
    if not pdb_path.exists():
        raise RuntimeError(f"Input file does not exist: `{pdb_path.absolute()}`")
    if not dcd_path.exists():
        raise RuntimeError(f"Input file does not exist: `{dcd_path.absolute()}`")

    # Set template notebook and user notebook filepath
    notebook_path = Path(_version.__file__).parent / "../docs/tutorials/dynophore.ipynb"
    new_notebook_path = Path(workspace_path) / "dynophore.ipynb"

    # Copy template notebook to user-defined workspace
    if notebook_path.exists():
        print("\nCopy dynophore notebook to user workspace...")
        copyfile(notebook_path, new_notebook_path)
        if new_notebook_path.exists():
            print(f"Dynophore notebook location: `{new_notebook_path.absolute()}`")
        else:
            raise RuntimeError(
                f"Could not create dynophore notebook at selected location f"
                f"`{new_notebook_path.absolute()}`"
            )
    else:
        raise RuntimeError(
            f"Could not find dynophore notebook at expected location `{notebook_path.absolute()}`."
        )

    # Replace template filepaths in notebook with user-defined filepaths
    print("\nUpdate filepaths in notebook to user filepaths...")
    search_replace_tuples = [
        ("../../dynophores/tests/data/1KE7-1/DynophoreApp", str(dyno_path.absolute())),
        ("../../dynophores/tests/data/1KE7-1/startframe.pdb", str(pdb_path.absolute())),
        ("../../dynophores/tests/data/1KE7-1/trajectory.dcd", str(dcd_path.absolute())),
    ]
    _update_paths_in_notebook(new_notebook_path, search_replace_tuples)


def _open_notebook(notebook_path):
    """
    Open input notebook with Jupyter Lab.

    Notes
    -----
    Same behaviour like `jupyter lab notebook_path`.
    """

    notebook_path = Path(notebook_path)
    if not notebook_path.exists():
        raise RuntimeError(f"Input file does not exist: `{notebook_path.absolute()}`")

    print("Open dynophore notebook with Jupyter Lab...")
    subprocess.run(["jupyter", "lab", notebook_path.absolute()])


def _update_paths_in_notebook(notebook_path, search_replace_tuples):
    """
    Read notebook file as string, replace all search instances (filepaths), and overwrite notebook
    file with this updated string.
    """

    # Read in the file
    with open(notebook_path, "r") as f:
        filedata = f.read()

    # Replace the target string
    for (search_str, replace_str) in search_replace_tuples:
        filedata = filedata.replace(search_str, replace_str)

    # Write the file out again
    with open(notebook_path, "w") as f:
        f.write(filedata)
