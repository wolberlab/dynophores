"""
Command Line Interface for the project.
"""

import argparse
from pathlib import Path
from shutil import copyfile
import subprocess

from . import _version
from . import api


PATH_TEST_DATA = Path(__file__).parent / "tests" / "data"


def main():
    """
    Main CLI function for the following commands:
    - dynophore create
    - dynophore visualize
    - dynophore demo
    """

    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers()

    create_subparser = subparsers.add_parser("create")
    visualize_subparser = subparsers.add_parser("visualize")
    demo_subparser = subparsers.add_parser("demo")

    # Arguments and function to balled for sub-command create
    create_subparser.add_argument(
        "-p",
        "--pmz",
        type=str,
        help="Path to input .PDB or .PMZ file containing a PDB or LigandScout binding site",
        required=True,
    )
    create_subparser.add_argument(
        "-d",
        "--dcd",
        type=str,
        help="Path to input .DCD file containing the molecular dynamics trajectory",
        required=True,
    )
    create_subparser.add_argument(
        "-o",
        "--out",
        type=str,
        help="Path to output location for DynophoreApp output folder, must end with /",
        required=True,
    )
    create_subparser.add_argument(
        "-n", "--name", type=str, help="Name for dynophore", required=True
    )
    create_subparser.add_argument(
        "-f",
        "--feature-def-file",
        type=str,
        help="Chemical feature definitions file ",
        required=False,
    )
    create_subparser.add_argument(
        "-3", "--three-letter-code", type=str, help="Ligand 3-letter code", required=False
    )
    create_subparser.add_argument("-c", "--chain", type=str, help="Ligand chain", required=False)
    create_subparser.set_defaults(func=_create)

    # Arguments and function to be called for sub-command visualize
    visualize_subparser.add_argument(
        "-i",
        "--dyno",
        type=str,
        help="Path to DynophoreApp output folder",
        required=True,
    )
    visualize_subparser.add_argument(
        "-p",
        "--pdb",
        type=str,
        help="Path to pdb (topology) file from trajectory, e.g. first frame",
        required=True,
    )
    visualize_subparser.add_argument(
        "-d",
        "--dcd",
        type=str,
        help="Path to dcd (trajectory) file",
        required=False,
    )
    visualize_subparser.set_defaults(func=_visualize)

    # Arguments and function to be called for sub-command demo
    demo_subparser.add_argument(
        "out",
        type=str,
        help="Path to target output folder",
    )
    demo_subparser.set_defaults(func=_demo)

    args = parser.parse_args()
    try:
        args.func(args)
    except AttributeError:
        # Run help if no arguments were given
        subprocess.run(["dynophore", "-h"])


def _greet():
    """
    Print CLI greeting.
    """

    logo_path = Path(_version.__file__).parent / "data/stegosaurus.txt"
    with open(logo_path, "r", encoding="ascii") as f:
        print(f.read())
    title_str = f"Dynophores CLI {_version.get_versions()['version']}"
    print(f"\n{title_str:^64}")


def _create(args):
    """
    Helper function for CLI to create dynophore data: `dynophore create` subcommand
    """

    pmz_or_pdb_path = Path(args.pmz)
    dcd_path = Path(args.dcd)
    out_path = Path(args.out)
    name = args.name
    if args.feature_def_file is not None:
        feature_def_path = Path(args.feature_def_file)
    else:
        feature_def_path = None
    three_letter_code = args.three_letter_code
    chain = args.chain

    api.create.create(
        pmz_or_pdb_path, dcd_path, out_path, name, feature_def_path, three_letter_code, chain
    )


def _visualize(args):
    """
    Helper function for CLI to visualize dynophore data: `dynophore visualize` subcommand
    """

    dyno_path = Path(args.dyno)
    pdb_path = Path(args.pdb)
    if args.dcd is not None:
        dcd_path = Path(args.dcd)
    else:
        dcd_path = None

    api.visualize.visualize(dyno_path, pdb_path, dcd_path)


def _demo(args):
    """
    Helper function for CLI to visualize dynophore data: `dynophore demo` subcommand
    """

    out_path = Path(args.out)
    api.demo.demo(out_path)
