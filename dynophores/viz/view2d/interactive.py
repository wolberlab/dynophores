"""
Contains RDKit ligand 2D visualizations (interactive).
"""

from ipywidgets import interact, fixed
import ipywidgets as widgets

from dynophores.viz import view2d


def show(dynophore):
    """
    Show interactive ligand 2D view.
    Toogle between superfeatures on/off and atom serial numbers on/off.

    Parameters
    ----------
    dynophore : dynophores.Dynophore
        Dynophore object.
    """

    interact(
        view2d.static.show,
        dynophore=fixed(dynophore),
        show_superfeatures=widgets.Checkbox(value=False, description="Show superfeatures"),
        show_pdb_serial_numbers=widgets.Checkbox(
            value=False, description="Show atom serial numbers"
        ),
    )
