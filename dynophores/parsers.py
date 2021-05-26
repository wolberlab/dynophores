"""
Contains parsers for
- Dynophore JSON file, which contains statistis on the occurrences and distances of superfeatures
  and environmental partners
- Dynophore PML file, which contains 3D coordinates for points in superfeature
  point clouds
"""

import json
from pathlib import Path
import xml.etree.ElementTree as ET

import numpy as np


def _json_to_dict(json_path):
    """
    Parse JSON file content.

    Parameters
    ----------
    json_path : str or pathlib.Path
        Path to JSON file.

    Returns
    -------
    dict
        Dynophore data with the following keys and nested keys (key : value data type):
        - id : str
        - ligand_name : str
        - ligand_smiles : str
        - ligand_mdl_mol_buffer : str
        - ligand_atom_serials : str
        - superfeatures : list
            - id : str
            - feature_type : str
            - atom_numbers : list of int
            - occurrences : list of int
            - envpartners : list
                - id : str
                - name : str
                - atom_numbers : list of int
                - occurrences : list of int
                - distances : list of float
    """

    with open(json_path, "r") as f:
        json_string = f.read()
        dynophore_dict = json.loads(json_string)

    return dynophore_dict


def _pml_to_dict(pml_path):
    """
    Parse PML file content (selections of it).

    Parameters
    ----------
    pml_path : str or pathlib.Path
        Path to PML file.

    Returns
    -------
    dict
        Superfeature data with the following keys and nested keys (key : value data type):

        Example:
        - <superfeature id>
          - id : str
          - color : str
          - center : numpy.array
          - points : list
            - x : float
            - y : float
            - z : float
            - frame_ix : int
            - weight : float
    """

    dynophore3d_xml = ET.parse(pml_path)
    dynophore3d_dict = {}

    feature_clouds = dynophore3d_xml.findall("featureCloud")
    for feature_cloud in feature_clouds:

        # Superfeature ID
        superfeature_feature_name = feature_cloud.get("name")
        superfeature_atom_numbers = feature_cloud.get("involvedAtomSerials")
        superfeature_id = f"{superfeature_feature_name}[{superfeature_atom_numbers}]"
        # Superfeature color
        superfeature_color = feature_cloud.get("featureColor")
        # Superfeature cloud center
        center = feature_cloud.find("position")
        center_data = np.array(
            [
                float(center.get("x3")),
                float(center.get("y3")),
                float(center.get("z3")),
            ]
        )
        # Superfeature cloud points
        additional_points = feature_cloud.findall("additionalPoint")
        additional_points_data = []
        for additional_point in additional_points:
            additional_point_data = {
                "x": float(additional_point.get("x3")),
                "y": float(additional_point.get("y3")),
                "z": float(additional_point.get("z3")),
                "frame_ix": int(additional_point.get("frameIndex")),
                "weight": float(additional_point.get("weight")),
            }
            additional_points_data.append(additional_point_data)

        dynophore3d_dict[superfeature_id] = {}
        dynophore3d_dict[superfeature_id]["id"] = superfeature_id
        dynophore3d_dict[superfeature_id]["color"] = superfeature_color
        dynophore3d_dict[superfeature_id]["center"] = center_data
        dynophore3d_dict[superfeature_id]["points"] = additional_points_data

    return dynophore3d_dict


def _json_pml_to_dict(json_path, pml_path):
    """
    Parse JSON and PML file content into one dictionary.

    Parameters
    ----------
    json_path : str or pathlib.Path
        Path to JSON file.
    pml_path : str or pathlib.Path
        Path to PML file.

    Returns
    -------
    dict
        Dynophore data with the following keys and nested keys (key : value data type):
        - "id" : str
        - "ligand" : dict
          - "name" : str
          - "smiles" : str
          - "mdl_mol_buffer" : str
          - "atoms_serials" : list of int
        - "superfeatures" : dict
            - <superfeature id> : str
                - "id" : str
                - "feature_type" : str
                - "atom_numbers" : list of int
                - "occurrences" : list of int
                - "envpartners" : dict
                    - <envpartner id> : str
                        - "id" : str
                        - "atom_numbers" : list of int
                        - "occurrences" : list of int
                        - "distances" : numpy.array
                        - "residue_name" : str
                        - "residue_number" : int
                        - "chain" : int
                - "color" : str
                - "cloud" : dict
                    - "center" : numpy.array
                    - "points" : list
                        - "x" : float
                        - "y" : float
                        - "z" : float
                        - "frame_ix" : int
                        - "weight" : float
    """

    json_path = Path(json_path)
    pml_path = Path(pml_path)

    json_dict = _json_to_dict(json_path)
    pml_dict = _pml_to_dict(pml_path)

    # Check if superfeatures are the same in both files
    json_superfeatures = sorted([sp["id"] for sp in json_dict["superfeatures"]])
    pml_superfeatures = sorted(list(pml_dict.keys()))
    if json_superfeatures != pml_superfeatures:
        raise ValueError(
            f"Your PML and JSON files are not matching. Superfeatures must be the same.\n"
            f"Your JSON file: {json_path}\n"
            f"Your PML file: {pml_path}"
        )

    # Merge dicts from JSON and PML files
    dynophore_dict = {}
    dynophore_dict["id"] = json_dict["id"]
    dynophore_dict["ligand"] = _format_ligand_dict(json_dict)
    dynophore_dict["superfeatures"] = _format_superfeatures_dict(json_dict, pml_dict)

    return dynophore_dict


def _format_ligand_dict(json_dict):
    """
    Format ligand dictionary.

    Parameters
    ----------
    json_dict : dict
        Dictionary with JSON file information.

    Returns
    -------
    dict
        Formatted dictionary as needed to initialize a Ligand object using **kwargs.
    """

    ligand_dict = {}
    ligand_dict["name"] = json_dict["ligand_name"]
    ligand_dict["smiles"] = json_dict["ligand_smiles"]
    ligand_dict["mdl_mol_buffer"] = json_dict["ligand_mdl_mol_buffer"]
    ligand_dict["atom_serials"] = [
        int(serial) for serial in json_dict["ligand_atom_serials"].split(",")
    ]

    return ligand_dict


def _format_superfeatures_dict(json_dict, pml_dict):
    """
    Format superfeatures dictionary.

    Parameters
    ----------
    json_dict : dict
        Dictionary with JSON file information.
    pml_dict : dict
        Dictionary with PML file information.

    Returns
    -------
    dict
        Formatted dictionary as needed to initialize the dynophore's SuperFeature object using
        **kwargs.
    """

    # Save superfeatures as dict with superfeature IDs as keys
    superfeatures_dict = {
        superfeature["id"]: superfeature for superfeature in json_dict["superfeatures"]
    }

    for superfeature_id, superfeature in superfeatures_dict.items():
        # Add color to superfeatures
        color = pml_dict[superfeature_id]["color"]
        superfeature["color"] = color
        # Add cloud to superfeatures
        pml_dict[superfeature_id].pop("id")
        pml_dict[superfeature_id].pop("color")
        cloud = pml_dict[superfeature_id]
        superfeature["cloud"] = cloud
        # Save envpartners as dict with envpartner IDs as keys
        for envpartner in superfeature["envpartners"]:
            envpartner["residue_name"] = envpartner["name"].split("_")[0]
            envpartner["residue_number"] = int(envpartner["name"].split("_")[1])
            envpartner["chain"] = envpartner["name"].split("_")[2]
            envpartner["id"] = envpartner["id"].replace("_", "-")
            envpartner.pop("name")
        superfeature["envpartners"] = {
            envpartner["id"]: envpartner for envpartner in superfeature["envpartners"]
        }

    return superfeatures_dict
