"""
dynophores.core.ligand

Handles the Ligand class, which contains superfeature-focused ligand data.
"""

from collections import defaultdict

from rdkit import Chem, Geometry
from rdkit.Chem import AllChem
from rdkit.Chem.Draw import rdMolDraw2D, IPythonConsole
from IPython.display import Image

from dynophores.utils import pdb_ligand_data_for_rdkit, hex_to_rgb

# Set width/height of molecule view
IMAGE_WIDTH = 700
IMAGE_HEIGHT = 400

IPythonConsole.molSize = IMAGE_WIDTH, IMAGE_HEIGHT


class Ligand:
    """
    Class to store superfeature-focused ligand data.

    Attributes
    ----------
    name : str
        3-letter name of structure-bound ligand.
    smiles : str
        SMILES of structure-bound ligand.
    _superfeatures_pdb_atom_ids : dict of str: list of int
        Per superfeature ID (keys): involved atom PDB IDs of ligand.
    _superfeatures_colors : dict of str: tuple of float
        Per superfeature ID (keys): superfeature color as scaled RGB values.
    _pdb_block : str
        Ligand PDB block
    _pdb_atom_ids : list of int
        Ligand PDB atom IDs.
    """

    def __init__(
        self,
        name,
        smiles,
        superfeatures_pdb_atom_ids,
        superfeatures_colors,
        pdb_block,
        pdb_atom_ids,
        **kwargs
    ):

        self.name = name
        self.smiles = smiles
        self._superfeatures_pdb_atom_ids = superfeatures_pdb_atom_ids
        self._superfeatures_colors = superfeatures_colors
        self._pdb_block = pdb_block
        self._pdb_atom_ids = pdb_atom_ids

    @classmethod
    def from_pdb(cls, pdb_path, dynophore):
        """
        Load ligand data from PDB file and dynophore object.  # TODO add ligand PDB content to PML

        Parameters
        ----------
        pdb_path : str or pathlib.Path
            Path to PDB file
        dynophore : dynophores.Dynophore
            Dynophore object.

        Returns
        -------
        dynophores.core.ligand.Ligand
            Ligand object.
        """

        pdb_block, pdb_atom_ids = pdb_ligand_data_for_rdkit(pdb_path, dynophore.ligand_name)
        return cls(
            name=dynophore.ligand_name,
            smiles=dynophore.ligand_smiles,
            superfeatures_pdb_atom_ids={
                superfeature_id: superfeature.atom_numbers
                for superfeature_id, superfeature in dynophore.superfeatures.items()
            },
            superfeatures_colors={
                superfeature_id: tuple(hex_to_rgb(superfeature.color))
                for superfeature_id, superfeature in dynophore.superfeatures.items()
            },
            pdb_block=pdb_block,
            pdb_atom_ids=pdb_atom_ids,
        )

    def view2d(self, show_pdb_serial_numbers=False):
        """
        2D view of ligand. Uses rdkit.

        Parameters
        ----------
        show_pdb_serial_numbers : bool
            Show PDB serial numbers of ligand atoms (default: False).

        Returns
        -------


        Notes
        -----
        For now, the procedure is a bit complicated to get from the ligand PDB block to the rdkit
        ligand. PDB atom serials must be set as properties to the rdkit molecule and bond orders
        must be inferred from SMILES.
        TODO Wait for MDAnalysis-RDKit converter:
        https://www.mdanalysis.org/2020/08/29/gsoc-report-cbouy/
        https://github.com/MDAnalysis/mdanalysis/issues/2468
        """

        # Read ligand from PDB block
        mol = Chem.MolFromPDBBlock(self._pdb_block, removeHs=False)

        # Set PDB atom serial numbers to rdkit molecule
        if len(self._pdb_atom_ids) == mol.GetNumAtoms():
            for atom, pdb_atom_id in zip(mol.GetAtoms(), self._pdb_atom_ids):
                atom.SetIntProp("pdb_atom_id", int(pdb_atom_id))
        else:
            raise ValueError("Number of atoms in RDKit and MDAnalysis are not the same.")

        # Remove hydrogens
        mol = Chem.RemoveHs(mol)

        # Assign bond orders based on SMILES
        reference_mol = Chem.MolFromSmiles(self.smiles)
        mol = AllChem.AssignBondOrdersFromTemplate(reference_mol, mol)
        AllChem.Compute2DCoords(mol)

        if show_pdb_serial_numbers:
            for atom in mol.GetAtoms():
                atom.SetAtomMapNum(atom.GetIntProp("pdb_atom_id"))

        return mol

    def view2d_superfeatures(self, show_pdb_serial_numbers=False):
        """
        2D view of ligand with highlighted superfeatures. Uses rdkit.

        Notes
        -----
        Code to highlight substructures is taken from Greg Landrum's blog post:
        http://rdkit.blogspot.com/2020/10/molecule-highlighting-and-r-group.html
        """

        mol = self.view2d(show_pdb_serial_numbers)
        superfeatures_pdb_atom_ids = self._superfeatures_pdb_atom_ids
        superfeatures_colors = self._superfeatures_colors

        (
            highlight_atoms,
            highlight_bonds,
            highlight_radius,
            highlight_linewidth,
            rings,
        ) = self._get_superfeatures_drawing_data(
            mol, superfeatures_pdb_atom_ids, superfeatures_colors
        )
        img = self._draw_superfeatures(
            mol, highlight_atoms, highlight_bonds, highlight_radius, highlight_linewidth, rings
        )

        return img

    @staticmethod
    def _get_superfeatures_drawing_data(mol, superfeatures_pdb_atom_ids, superfeatures_colors):
        """
        Get atoms/bonds/radius/linewidth/radius needed to highlight superfeatures in ligand.

        Parameters
        ----------
        mol : rdkit.Chem.rdchem.Mol
            Molecule.
        superfeatures_pdb_atom_ids : dict of str: list of int
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

        for superfeature_id, pdb_atom_ids in superfeatures_pdb_atom_ids.items():

            # Map PDB atom IDs to RDKit atom IDs
            rdkit_atom_idx = [
                atom.GetIdx()
                for atom in mol.GetAtoms()
                if atom.GetIntProp("pdb_atom_id") in pdb_atom_ids
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

    @staticmethod
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
        for (ring, color) in rings:
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
