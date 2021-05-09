"""
Contains NGLview 3D visualizations.
"""

import warnings

import nglview as nv
import MDAnalysis as mda

from dynophores.utils import hex_to_rgb

warnings.filterwarnings("ignore", category=DeprecationWarning)

MACROMOLECULE_COLOR = "#005780"  # blue
LIGAND_COLOR = "#808080"  # grey


def show(dynophore, pdb_path, dcd_path=None, visualization_type="spheres", macromolecule_color=MACROMOLECULE_COLOR, ligand_color=LIGAND_COLOR):
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
    envpartners = " or ".join([f"(:{chain} and {residue_number})" for chain, residue_number in envpartners])
    view.add_representation("hyperball", selection=envpartners)

    # Add dynophore
    _add_dynophore(view, dynophore, visualization_type)

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


def _add_dynophore(view, dynophore, visualization_type):
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
        buffer = {"position": [], "color": [], "radius": []}
        for point in superfeature.cloud.points:
            buffer["position"] += [point.x, point.y, point.z]
            buffer["color"] += hex_to_rgb(superfeature.color)
            buffer["radius"] += [0.1]
        
        if visualization_type == "spheres":
            js = _js_sphere_buffer(buffer, superfeature.id)
        elif visualization_type == "points":
            js = _js_point_buffer(buffer, superfeature.id)
        else:
            raise ValueError("Visualization style unknown. Please choose from `spheres` or `points`.")

        view._js(js)

    return view

def _js_sphere_buffer(buffer, superfeature_id):

    return f"""
        var params = {buffer};
        var shape = new NGL.Shape('{superfeature_id}');
        var buffer = new NGL.SphereBuffer(params);

        shape.addBuffer(buffer);
        var shapeComp = this.stage.addComponentFromObject(shape);
        shapeComp.addRepresentation("buffer", {{opacity:0.55}});
        """

def _js_point_buffer(buffer, superfeature_id):

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
        shapeComp.addRepresentation("buffer", {{ opacity: 0.55 }});
        """
