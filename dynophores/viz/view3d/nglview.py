"""
Contains NGLview visualizations.
"""

import warnings

import nglview as nv
import MDAnalysis as mda
import matplotlib

from ...parsers import _pml_to_dict
from ...definitions import FEATURE_COLORS

warnings.filterwarnings("ignore", category=DeprecationWarning)


def show_dynophore3d(pml_path, pdb_path, dcd_path=None):

    if dcd_path is None:
        view = _show_structure(pdb_path)
    else:
        view = _show_trajectory(pdb_path, dcd_path)

    view.clear_representations()
    view.add_representation("cartoon", selection="protein", color="grey")
    view.add_representation("ball+stick", selection="ligand")

    view = _add_dynophore(view, pml_path)
    return view


def _show_structure(pdb_path):

    view = nv.show_file(str(pdb_path), ext="pdb")
    return view


def _show_trajectory(pdb_path, dcd_path):

    md_universe = mda.Universe(str(pdb_path), str(dcd_path))
    view = nv.show_mdanalysis(md_universe)
    return view


def _add_dynophore(view, pml_path):

    dynophore3d_dict = _pml_to_dict(pml_path)

    for superfeature_id, cloud in dynophore3d_dict.items():
        sphere_buffer = {"position": [], "color": [], "radius": []}
        # TODO: Add cloud name!
        for point_coordinates in cloud["coordinates"]:
            sphere_buffer["position"] += point_coordinates.tolist()
            sphere_buffer["color"] += matplotlib.colors.to_rgb(
                FEATURE_COLORS[superfeature_id.split("[")[0]]
            )
            sphere_buffer["radius"] += [0.1]
        view.shape.add_buffer("sphere", **sphere_buffer)

    return view
