"""
Contains parser for LigandScout pml file, which contains 3D coordinates for points in superfeature
point clouds.
"""

import xml.etree.ElementTree as ET


def _parse_pml(filepath):
    """
    Parse pml file content.

    Parameters
    ----------
    filepath : str or pathlib.Path
        Path to pml file.

    Returns
    -------
    dict of (str: list of list of float)
        Per superfeature, coordinates for points in point cloud.
    """

    dynophore3d_xml = ET.parse(filepath)
    dynophore3d_dict = {}

    feature_clouds = dynophore3d_xml.findall("featureCloud")
    for feature_cloud in feature_clouds:

        superfeature_id = f"{feature_cloud.get('name')}-{feature_cloud.get('featureId')}"
        dynophore3d_dict[superfeature_id] = []

        position = feature_cloud.find("position")
        position_coordinates = [i[1] for i in position.items()[:3]]
        dynophore3d_dict[superfeature_id].append(position_coordinates)

        additional_points = feature_cloud.findall("additionalPoint")
        for additional_point in additional_points:
            additional_point_coordinates = [i[1] for i in additional_point.items()[:3]]
            dynophore3d_dict[superfeature_id].append(additional_point_coordinates)

    return dynophore3d_dict
