"""
Contains NGLview 3D visualizations.
"""

import warnings

import numpy as np
import nglview as nv
import MDAnalysis as mda
from matplotlib import colors

from dynophores.utils import hex_to_rgb_saturation_sequence

warnings.filterwarnings("ignore", category=DeprecationWarning)

MACROMOLECULE_COLOR = "#005780"  # blue
LIGAND_COLOR = "#808080"  # grey


def show(
    dynophore,
    pdb_path,
    dcd_path=None,
    visualization_type="spheres",
    color_cloud_by_frame=False,
    macromolecule_color=MACROMOLECULE_COLOR,
):
    """
    Show the dynophore point cloud with its ligand-bound structure and optionally the underlying
    MD trajectory.

    Parameters
    ----------
    dynophore : dynophore.Dynophore
        Dynophore.
    pdb_path : str or pathlib.Path
        Path to PDB file (structure/topology)
    dcd_path : None or str or pathlib.Path
        Optionally: Path to DCD file (trajectory).
    visualization_type : str
        Visualization type for dynophore cloud: `spheres` or `points`
    color_cloud_by_frame : bool
        Use color saturation to indicate frame index of each cloud point (default: False).
    macromolecule_color : str
        Hex code for macromolecule color.

    Returns
    -------
    nglview.widget.NGLWidget
        Visualization with the NGL Viewer.
    """

    # Show structure or trajectory
    if dcd_path is None:
        view = _show_structure(pdb_path)
    else:
        view = _show_trajectory(pdb_path, dcd_path)
    # Set representation and color for protein and ligand
    view.clear_representations()
    view.add_representation("cartoon", selection="protein", color=macromolecule_color)
    view.add_representation("hyperball", selection="ligand")

    # Add interacting pocket residues
    envpartners = dynophore.unique_envpartners_chain_residue_number
    envpartners_string = " or ".join(
        [f"(:{chain} and {residue_number})" for chain, residue_number in envpartners]
    )
    view.add_representation("licorice", selection=envpartners_string)
    # Add dynophore
    _add_dynophore(view, dynophore, visualization_type, color_cloud_by_frame)

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
    pdb_path : str or pathlib.Path
        Path to PDB file (structure/topology)
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


def _add_dynophore(view, dynophore, visualization_type, color_cloud_by_frame=False):
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
    visualization_type : str
        Visualization type for dynophore cloud: `spheres` or `points`
    color_cloud_by_frame : bool
        Use color saturation to indicate frame index of each cloud point (default: False).


    Returns
    -------
    nglview.widget.NGLWidget
        Visualization with the NGL Viewer.
    """

    for _, superfeature in dynophore.superfeatures.items():
        buffer = {"position": [], "color": [], "radius": []}
        cloud = superfeature.cloud.data.copy()

        # Use color saturation to indicate frame index per point
        if color_cloud_by_frame:
            colors_by_frame = hex_to_rgb_saturation_sequence(
                f"#{superfeature.color}", dynophore.n_frames
            )
            color_matrix = colors_by_frame[cloud["frame_ix"].to_list()]
            cloud["rgb"] = color_matrix.tolist()
            opacity = 0.8
        # Use the same color for all points of the same feature type
        else:
            color = colors.hex2color(f"#{superfeature.color}")
            color_matrix = np.repeat(np.array([color]), [len(cloud)], axis=0)
            cloud["rgb"] = color_matrix.tolist()
            opacity = 0.55

        for _, point in cloud.iterrows():
            buffer["position"] += [point["x"], point["y"], point["z"]]
            buffer["radius"] += [0.1]
            buffer["color"] += point["rgb"]

        if visualization_type == "spheres":
            js = _js_sphere_buffer(buffer, superfeature.id, opacity)
        elif visualization_type == "points":
            js = _js_point_buffer(buffer, superfeature.id, opacity)
        else:
            raise ValueError(
                "Visualization style unknown. Please choose from `spheres` or `points`."
            )

        view._js(js)

    return view


def _js_sphere_buffer(buffer, superfeature_id, opacity):
    """
    Create JavaScript string generating a sphere buffer from buffer parameters.

    Parameters
    ----------
    buffer : str
        Buffer parameters.
    superfeature_id : str
        Superfeature ID.
    opacity : float
        Sphere opacity (the higher the less transparent).

    Returns
    -------
    JavaScript string generating a sphere buffer.
    """

    return f"""
        var params = {buffer};
        var shape = new NGL.Shape('{superfeature_id}');
        var buffer = new NGL.SphereBuffer(params);

        shape.addBuffer(buffer);
        var shapeComp = this.stage.addComponentFromObject(shape);
        shapeComp.addRepresentation("buffer", {{opacity:{opacity}}});
        """


def _js_point_buffer(buffer, superfeature_id, opacity):
    """
    Create JavaScript string generating a point buffer from buffer parameters.
    TODO: Visualization contains artifacts:
    - https://github.com/dominiquesydow/dynophores/issues/27
    - https://github.com/nglviewer/ngl/issues/868

    Parameters
    ----------
    buffer : str
        Buffer parameters.
    superfeature_id : str
        Superfeature ID.
    opacity : float
        Sphere opacity (the higher the less transparent).

    Returns
    -------
    JavaScript string generating a point buffer.
    """

    return f"""
        var point_buffer = new NGL.PointBuffer(
            {{
                position: new Float32Array({buffer["position"]}),
                color: new Float32Array({buffer["color"]})
            }},
            {{
                sizeAttenuation: true,
                pointSize: 2,
                opacity: 0.1,
                useTexture: true,
                alphaTest: 0.0,
                edgeBleach: 0.7,
                forceTransparent: true,
                sortParticles: true
            }}
        )
        var shape = new NGL.Shape('{superfeature_id}');

        shape.addBuffer(point_buffer);
        var shapeComp = this.stage.addComponentFromObject(shape);
        shapeComp.addRepresentation("buffer", {{ opacity: {opacity} }});
        """
