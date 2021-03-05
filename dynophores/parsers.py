"""
Contains parsers for
- Dynophore JSON file, which contains statistis on the occurrences and distances of superfeatures
  and environmental partners
- Dynophore PML file, which contains 3D coordinates for points in superfeature
  point clouds
"""

import json
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
        TODO
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
    dict of (str: list of list of float)
        Per superfeature, coordinates for points in point cloud.
    """

    dynophore3d_xml = ET.parse(pml_path)
    dynophore3d_dict = {}

    feature_clouds = dynophore3d_xml.findall("featureCloud")
    for feature_cloud in feature_clouds:

        superfeature_feature_name = feature_cloud.get("name")
        superfeature_atom_numbers = feature_cloud.get("involvedAtomSerials")
        superfeature_id = f"{superfeature_feature_name}[{superfeature_atom_numbers}]"

        dynophore3d_dict[superfeature_id] = {}
        dynophore3d_dict[superfeature_id]["id"] = superfeature_id

        center = feature_cloud.find("position")
        center_coordinates = [float(i[1]) for i in center.items()[:3]]
        dynophore3d_dict[superfeature_id]["center"] = np.array(center_coordinates)

        additional_points = feature_cloud.findall("additionalPoint")
        additional_point_coordinates = np.array(
            [
                [float(i[1]) for i in additional_point.items()[:3]]
                for additional_point in additional_points
            ]
        )

        dynophore3d_dict[superfeature_id]["points"] = additional_point_coordinates

    return dynophore3d_dict
