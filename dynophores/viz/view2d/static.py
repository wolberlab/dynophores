"""
Contains RDKit ligand 2D visualizations (static).
"""
from collections import defaultdict

from rdkit import Geometry
from rdkit.Chem.Draw import rdMolDraw2D, IPythonConsole
from IPython.display import Image
from matplotlib import colors

# Set width/height of molecule view
IMAGE_WIDTH = 700
IMAGE_HEIGHT = 400

IPythonConsole.molSize = IMAGE_WIDTH, IMAGE_HEIGHT


def show(dynophore, show_superfeatures=False, show_pdb_serial_numbers=False):
    """
    Show static ligand 2D view.

    Parameters
    ----------
    dynophore : dynophores.Dynophore
        Dynophore object.

    Returns
    -------
    rdkit.Chem.rdchem.Mol or IPython.core.display.Image
        Molecule (with/without superfeatures).
    """

    mol = dynophore.ligand.rdkit_molecule(show_pdb_serial_numbers)

    if show_superfeatures:
        return _add_superfeatures(mol, dynophore)
    else:
        return mol


def _add_superfeatures(mol, dynophore):
    """
    2D view of ligand with highlighted superfeatures. Uses rdkit.

    Parameters
    ----------
    mol : rdkit.Chem.rdchem.Mol or IPython.core.display.Image
        Molecule (without superfeatures).
    dynophore : dynophores.Dynophore
        Dynophore object.

    Returns
    -------
    IPython.core.display.Image
        2D view of ligand with superfeatures.

    Notes
    -----
    Code to highlight substructures is taken from Greg Landrum's blog post:
    http://rdkit.blogspot.com/2020/10/molecule-highlighting-and-r-group.html
    """

    superfeatures_atom_serials = dynophore.superfeatures_atom_serials
    superfeatures_colors = {
        superfeature_id: tuple(colors.hex2color(f"#{color}"))
        for superfeature_id, color in dynophore.superfeatures_colors.items()
    }

    (
        highlight_atoms,
        highlight_bonds,
        highlight_radius,
        highlight_linewidth,
        rings,
    ) = _get_superfeatures_drawing_data(mol, superfeatures_atom_serials, superfeatures_colors)
    img = _draw_superfeatures(
        mol, highlight_atoms, highlight_bonds, highlight_radius, highlight_linewidth, rings
    )

    return img


def _get_superfeatures_drawing_data(mol, superfeatures_atom_serials, superfeatures_colors):
    """
    Get atoms/bonds/radius/linewidth/radius needed to highlight superfeatures in ligand.

    Parameters
    ----------
    mol : rdkit.Chem.rdchem.Mol
        Molecule.
    superfeatures_atom_serials : dict of str: list of int
        Per superfeature ID (keys): involved atom PDB IDs of ligand.
    superfeatures_colors : dict of str: tuple of float
        Per superfeature ID (keys): superfeature color as scaled RGB values.

    Returns
    -------
    highlight_atoms : collections.defaultdict
        Atoms to be highlighted.
        Example:
        defaultdict(list, {22: [(0.969, 0.243, 0.243)], 7: [(0.18, 0.192, 0.573)]})
    highlight_bonds : collections.defaultdict
        Bonds to be highlighted.
        Example:
        defaultdict(list, {7: [(0.18, 0.192, 0.573), (0.18, 0.192, 0.573)]})
    highlight_radius : dict
        Atom highlight radius.
        Example:
        {22: 0.4, 7: 0.4}
    highlight_linewidth :
        Bond hightlight line width.
        Example:
        {7: 2, 11: 2}
    rings : tuple of tuple
        Rings to be filled/highlighted.
        Example:
        [
            ((8, 7, 11, 10, 9), (0.18, 0.192, 0.573)),
            ((18, 19, 25, 26, 27, 17), (0.18, 0.192, 0.573))
        ]

    Notes
    -----
    Making this a static method, so that it is clearer what the required input for this
    method it - but could use attributes instead.
    """

    highlight_atoms = defaultdict(list)
    highlight_bonds = defaultdict(list)
    highlight_radius = {}
    highlight_linewidth = {}
    rings = []

    mol_rings = mol.GetRingInfo().AtomRings()

    for superfeature_id, atom_serials in superfeatures_atom_serials.items():
        # Map PDB atom IDs to RDKit atom IDs
        rdkit_atom_idx = [
            atom.GetIdx()
            for atom in mol.GetAtoms()
            if atom.GetIntProp("pdb_atom_serial") in atom_serials
        ]
        # Get subpocket color
        color = superfeatures_colors[superfeature_id]

        # Highlight information for rings
        for ring in mol_rings:
            # Add only those rings that belong to current fragment
            if len(set(ring) - set(rdkit_atom_idx)) == 0:
                rings.append((ring, color))

        for atom_idx in rdkit_atom_idx:
            # Highlight information for atoms
            highlight_atoms[atom_idx].append(color)
            highlight_radius[atom_idx] = 0.4

            # Highlight information for bonds
            atom = mol.GetAtomWithIdx(atom_idx)
            bonds = atom.GetBonds()
            for bond in bonds:
                atom_idx_begin = bond.GetBeginAtom().GetIdx()
                atom_idx_end = bond.GetEndAtom().GetIdx()
                if (atom_idx_begin in rdkit_atom_idx) and (atom_idx_end in rdkit_atom_idx):
                    bond_idx = bond.GetIdx()
                    highlight_bonds[bond_idx].append(color)
                    highlight_linewidth[bond_idx] = 2

    return highlight_atoms, highlight_bonds, highlight_radius, highlight_linewidth, rings


def _draw_superfeatures(
    mol, highlight_atoms, highlight_bonds, highlight_radius, highlight_linewidth, rings
):
    """
    Draw superfeatures on 2D view of ligand.
    Use only with _get_superfeatures_drawing_data().

    Parameters
    ----------
    mol : rdkit.Chem.rdchem.Mol
        Molecule.
    highlight_atoms : collections.defaultdict
        Atoms to be highlighted.
    highlight_bonds : collections.defaultdict
        Bonds to be highlighted.
    highlight_radius : dict
        Atom highlight radius.
    highlight_linewidth :
        Bond hightlight line width.
    rings : tuple of tuple
        Rings to be filled/highlighted.

    Returns
    -------
    IPython.core.display.Image
        2D view of ligand with superfeatures.
    """

    legend = ""
    width = IMAGE_WIDTH
    height = IMAGE_HEIGHT
    d2d = rdMolDraw2D.MolDraw2DCairo(width, height)
    dos = d2d.drawOptions()
    dos.useBWAtomPalette()

    # First draw circles to highlight full rings
    d2d.DrawMoleculeWithHighlights(
        mol,
        legend,
        dict(highlight_atoms),
        dict(highlight_bonds),
        highlight_radius,
        highlight_linewidth,
    )
    d2d.ClearDrawing()
    conf = mol.GetConformer()
    for ring, color in rings:
        positions = []
        for atom_idx in ring:
            position = Geometry.Point2D(conf.GetAtomPosition(atom_idx))
            positions.append(position)
        d2d.SetFillPolys(True)
        d2d.SetColour(color)
        d2d.DrawPolygon(positions)
    dos.clearBackground = False

    # Now draw on top (again) the molecule with highlighted atoms and bonds
    d2d.DrawMoleculeWithHighlights(
        mol,
        legend,
        dict(highlight_atoms),
        dict(highlight_bonds),
        highlight_radius,
        highlight_linewidth,
    )
    d2d.FinishDrawing()

    # Show PNG text as image
    img = Image(d2d.GetDrawingText())

    return img
