"""
Contains NGLview visualizations.
"""

import warnings

import nglview as nv
import MDAnalysis as mda

from dynophores.utils import hex_to_rgb

warnings.filterwarnings("ignore", category=DeprecationWarning)


def show_dynophore3d(pml_path, pdb_path, dcd_path=None):
    """
    Show the dynophore point cloud with its ligand-bound structure and optionally the underlying
    MD trajectory.

    Parameters
    ----------
    pml_path : str or pathlib.Path
        Path to PML file (dynophore).
    pdb_path : str or pathlib.Path
        Path to PDB file (structure/topology)
    dcd_path : None or str or pathlib.Path
        Optionally: Path to DCD file (trajectory).

    Returns
    -------
    nglview.widget.NGLWidget
        Visualization with the NGL Viewer.
    """

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
    """
    Show the ligand-bound structure that belongs to the dynophore. Could be the starting frame of
    the MD simulation.

    Parameters
    ----------
    pdb_path : str or pathlib.Path
        Path to PDB file (structure/topology)

    Returns
    -------
    nglview.widget.NGLWidget
        Visualization with the NGL Viewer.
    """

    view = nv.show_file(str(pdb_path), ext="pdb")
    return view


def _show_trajectory(pdb_path, dcd_path):
    """
    Show the MD trajectory that is the bases for the dynophore.

    Parameters
    ----------
    dcd_path : None or str or pathlib.Path
        Optionally: Path to DCD file (trajectory).

    Returns
    -------
    nglview.widget.NGLWidget
        Visualization with the NGL Viewer.
    """

    md_universe = mda.Universe(str(pdb_path), str(dcd_path))
    view = nv.show_mdanalysis(md_universe)
    return view


def _add_dynophore(view, dynophore):
    """
    Add the dynophore point cloud to an existing view of its underlying structure (and optionally
    its trajectory).

    Parameters
    ----------
    view : nglview.widget.NGLWidget
        NGL Viewer containing the ligand-bound structure (optionally including the trajectory)
        belonging to the dynophore.
    dynophore : Dynophore
        Dynophore data (includes data from JSON and PML file).

    Returns
    -------
    nglview.widget.NGLWidget
        Visualization with the NGL Viewer.
    """

    for _, superfeature in dynophore.superfeatures.items():
        sphere_buffer = {"position": [], "color": [], "radius": []}
        for point in superfeature.cloud.points:
            sphere_buffer["position"] += [point.x, point.y, point.z]
            sphere_buffer["color"] += hex_to_rgb(superfeature.color)
            sphere_buffer["radius"] += [0.1]
        js = f"""
        var params = {sphere_buffer};
        var shape = new NGL.Shape('{superfeature.id}');
        var buffer = new NGL.SphereBuffer(params);
        shape.addBuffer(buffer);
        var shapeComp = this.stage.addComponentFromObject(shape);
        shapeComp.addRepresentation("buffer");
        """
        view._js(js)

    return view
