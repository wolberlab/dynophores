"""
Command Line Interface for the project.
"""

import argparse
from pathlib import Path
from shutil import copyfile
import subprocess

from . import _version

PATH_TEST_DATA = Path(__name__).parent / "dynophores" / "dynophores" / "tests" / "data"


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

    dynoviz demo
    """

    _greet()

    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers()

    create_subparser = subparsers.add_parser("create")
    open_subparser = subparsers.add_parser("open")
    demo_subparser = subparsers.add_parser("demo")

    # Arguments and function to be called for sub-command create
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

    # Arguments and function to be called for sub-command open
    open_subparser.add_argument(
        "notebook",
        type=str,
        help="Path to dynophore notebook file",
    )
    open_subparser.set_defaults(func=_open_viz)

    # Arguments and function to be called for sub-command demo
    demo_subparser.add_argument(
        "workspace",
        type=str,
        help="Path to workspace folder",
    )
    demo_subparser.set_defaults(func=_demo_viz)

    args = parser.parse_args()
    args.func(args)


def _greet():
    """
    Print CLI greeting.
    """

    logo_path = Path(_version.__file__).parent / "data/stegosaurus.txt"
    with open(logo_path, "r", encoding="ascii") as f:
        print(f.read())
    title_str = f"Dynophores CLI {_version.get_versions()['version']}"
    print(f"\n{title_str:^64}")


def _create_viz(args):
    """
    Create visualization notebook based on command line arguments.
    """

    new_notebook_path = Path(args.workspace) / "dynophore.ipynb"
    _copy_notebook(new_notebook_path)
    _update_paths_in_notebook(new_notebook_path, args.dyno_path, args.pdb_path, args.dcd_path)
    _open_notebook(new_notebook_path)


def _open_viz(args):
    """
    Open visualization notebook based on command line arguments.
    """
    print("_open_viz")

    _open_notebook(args.notebook)


def _demo_viz(args):
    """
    Create and open demo visualization notebook based.
    """

    new_notebook_path = Path(args.workspace) / "dynophore_demo.ipynb"
    _copy_notebook(new_notebook_path)
    _update_paths_in_notebook(
        new_notebook_path,
        (PATH_TEST_DATA / "1KE7-1/DynophoreApp").absolute(),
        (PATH_TEST_DATA / "1KE7-1/startframe.pdb").absolute(),
        (PATH_TEST_DATA / "1KE7-1/trajectory.dcd").absolute(),
    )
    _open_notebook(new_notebook_path)


def _copy_notebook(new_notebook_path):
    """
    Copy template dynophore notebook to user-defined workspace.
    """

    # Set template notebook
    template_notebook_path = Path(_version.__file__).parent / "notebooks/dynophore.ipynb"

    # Set user notebook filepath
    new_notebook_path = Path(new_notebook_path)
    if new_notebook_path.suffix != ".ipynb":
        raise RuntimeError(
            f"Input file path must have suffix `.ipynb`. "
            f"Your input: `{new_notebook_path.absolute()}`"
        )
    if not new_notebook_path.exists():
        new_notebook_path.parent.mkdir(parents=True, exist_ok=True)

    # Copy template notebook to user-defined workspace
    if template_notebook_path.exists():
        print("\nCopy dynophore notebook to user workspace...")
        copyfile(template_notebook_path, new_notebook_path)
        if new_notebook_path.exists():
            print(f"Dynophore notebook location: `{new_notebook_path.absolute()}`")
        else:
            raise RuntimeError(
                f"Could not create dynophore notebook at selected location "
                f"`{new_notebook_path.absolute()}`"
            )
    else:
        raise RuntimeError(
            f"Could not find dynophore notebook at expected location "
            f"`{template_notebook_path.absolute()}`."
        )


def _update_paths_in_notebook(notebook_path, dyno_path, pdb_path, dcd_path):
    """
    Read notebook file as string, replace all search instances (filepaths), and overwrite notebook
    file with this updated string.
    """

    notebook_path = Path(notebook_path)
    dyno_path = Path(dyno_path)
    pdb_path = Path(pdb_path)
    dcd_path = Path(dcd_path)

    if not notebook_path.is_file():
        raise RuntimeError(f"Input is no file or does not exist: `{notebook_path.absolute()}`")
    if not dyno_path.is_dir():
        raise RuntimeError(f"Input is no file or file does not exist: `{dyno_path.absolute()}`")
    if not pdb_path.is_file():
        raise RuntimeError(f"Input is no file or file does not exist: `{pdb_path.absolute()}`")
    if not dcd_path.is_file():
        raise RuntimeError(f"Input is no file or file does not exist: `{dcd_path.absolute()}`")

    # Replace template filepaths in notebook with user-defined filepaths
    print("\nUpdate filepaths in notebook to user filepaths...")
    search_replace_tuples = [
        ("../tests/data/1KE7-1/DynophoreApp", str(dyno_path.absolute())),
        ("../tests/data/1KE7-1/startframe.pdb", str(pdb_path.absolute())),
        ("../tests/data/1KE7-1/trajectory.dcd", str(dcd_path.absolute())),
    ]

    # Read in the file
    with open(notebook_path, "r") as f:
        filedata = f.read()

    # Replace the target string
    for (search_str, replace_str) in search_replace_tuples:
        filedata = filedata.replace(search_str, replace_str)

    # Write the file out again
    with open(notebook_path, "w") as f:
        f.write(filedata)


def _open_notebook(notebook_path):
    """
    Open input notebook with Jupyter Lab.

    Notes
    -----
    Same behaviour like `jupyter lab notebook_path`.
    """

    print("_open_notebook")

    notebook_path = Path(notebook_path)
    print(notebook_path)

    if not notebook_path.exists():
        raise RuntimeError(f"Input path does not exist: `{notebook_path.absolute()}`")
    if not notebook_path.is_file():
        raise RuntimeError(f"Input path is not a file: `{notebook_path.absolute()}`")
    if notebook_path.suffix != ".ipynb":
        raise RuntimeError(
            f"Input file path must have suffix `.ipynb`. Your input: `{notebook_path.absolute()}`"
        )

    print("Open dynophore notebook with Jupyter Lab...")
    subprocess.run(["jupyter", "lab", str(notebook_path.absolute())])
