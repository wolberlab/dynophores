"""
Utility functions used to generate and visualize dynophore data.
"""

import json
from pathlib import Path
from shutil import copyfile


from .. import _version


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


def _update_paths_in_notebook(notebook_path, dyno_path, pdb_path, dcd_path=None):
    """
    Read notebook file as string, replace all search instances (filepaths), and overwrite notebook
    file with this updated string.
    """

    notebook_path = Path(notebook_path)
    if not notebook_path.is_file():
        raise RuntimeError(f"Input is no file or does not exist: `{notebook_path.absolute()}`")

    dyno_path = Path(dyno_path)
    if not dyno_path.is_dir():
        raise RuntimeError(f"Input is no file or file does not exist: `{dyno_path.absolute()}`")

    pdb_path = Path(pdb_path)
    if not pdb_path.is_file():
        raise RuntimeError(f"Input is no file or file does not exist: `{pdb_path.absolute()}`")

    if dcd_path is not None:
        dcd_path = Path(dcd_path)
        if not dcd_path.is_file():
            raise RuntimeError(f"Input is no file or file does not exist: `{dcd_path.absolute()}`")

    # Replace template filepaths in notebook with user-defined filepaths
    print("Update filepaths in notebook to user filepaths...")

    search_replace_tuples = [
        ("../tests/data/out", str(dyno_path.absolute())),
        ("../tests/data/in/startframe.pdb", str(pdb_path.absolute())),
    ]

    if dcd_path is not None:
        search_replace_tuples.append(("../tests/data/in/trajectory.dcd", str(dcd_path.absolute())))
    else:
        search_replace_tuples.append(
            (
                json.dumps('dcd_path = Path("../tests/data/in/trajectory.dcd")'),
                json.dumps("dcd_path = None"),
            )
        )

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
    Print message on how to open notebook with Jupyter Lab.
    """

    notebook_path = Path(notebook_path)

    print(
        f"""
\nFinished! To visualize the dynophore results run your notebook in Jupyter Lab:

    jupyter lab {notebook_path.absolute()}
    """
    )
