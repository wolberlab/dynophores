"""
Contains parsers for
- Dynophore JSON file, which contains statistis on the occurrences and distances of superfeatures
  and environmental partners
- Dynophore PML file, which contains 3D coordinates for points in superfeature
  point clouds
"""

from collections import OrderedDict

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
        Dynophore data with the following keys and nested keys (key : value data type):
        - id : str
        - ligand_name : str
        - ligand_smiles : str
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
    collections.OrderedDict of collections.OrderedDict
        Per superfeature cloud (key): ID, color, cloud center coordinates and point coordinates.

        Example:
        OrderedDict(
            [
                (
                    "H[4615,4623,4622,4613,4621,4614]",
                    OrderedDict(
                        [
                            ("id", "H[4615,4623,4622,4613,4621,4614]"),
                            ("color", "ffc20e"),
                            ("center", array([-14.740809, -6.0303836, -0.2289748])),
                            (
                                "points",
                                array(
                                    [
                                        [-15.363355, -3.327491, -2.418513],
                                        [-15.363355, -3.327491, -2.418513],
                                        [-14.156935, -7.772992, 0.909169],
                                    ]
                                ),
                            ),
                        ]
                    ),
                )
            ],
            [...]
        )
    """

    dynophore3d_xml = ET.parse(pml_path)
    dynophore3d_dict = OrderedDict()

    feature_clouds = dynophore3d_xml.findall("featureCloud")
    for feature_cloud in feature_clouds:

        superfeature_feature_name = feature_cloud.get("name")
        superfeature_atom_numbers = feature_cloud.get("involvedAtomSerials")
        superfeature_color = feature_cloud.get("featureColor")
        superfeature_id = f"{superfeature_feature_name}[{superfeature_atom_numbers}]"

        dynophore3d_dict[superfeature_id] = OrderedDict()
        dynophore3d_dict[superfeature_id]["id"] = superfeature_id
        dynophore3d_dict[superfeature_id]["color"] = superfeature_color

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
