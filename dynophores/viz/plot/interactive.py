"""
Contains interactive plotting functions for Jupyter notebooks.

IPyWidgets are amazing!

Resources:
- Great IPyWidgets documentation
  - https://ipywidgets.readthedocs.io/en/latest/examples/Using%20Interact.html#interactive
  - https://ipywidgets.readthedocs.io/en/latest/examples/Widget%20List.html
  - https://ipywidgets.readthedocs.io/en/latest/examples/Widget%20Styling.html
- Big thanks to @jaimergp for this getting-started notebook
  - https://github.com/volkamerlab/teachopencadd/blob/1fd1de46b66c5acbb4cdb8ef6083171cd84fb50f/
    teachopencadd/talktorials/T017_advanced_nglview_usage/talktorial.ipynb
"""

from ipywidgets import interact, fixed
import ipywidgets as widgets

from dynophores.viz import plot


def superfeatures_vs_envpartners(dynophore):
    """
    Plot heatmap of interactions between superfeatures and interaction partners.

    Parameters
    ----------
    dynophore : dynophores.Dynophore
        Dynophore.

    Returns
    -------
    function
        Parameterized IPyWidgets interact function.
    """

    style = {"description_width": "initial"}
    superfeature_ids = ["all"] + [superfeature.id for superfeature in dynophore.superfeatures]

    func = interact(
        plot.static.superfeatures_vs_envpartners,
        dynophore=fixed(dynophore),
        superfeature_names=widgets.SelectMultiple(
            options=superfeature_ids,
            value=["all"],
            description="Superfeature name(s):",
            style=style,
        ),
    )

    return func


def superfeatures_occurrences(dynophore):
    """
    Generate interactive widget to plot the superfeatures' occurrences as barcode.

    Parameters
    ----------
    dynophore : dynophores.Dynophore
        Dynophore.

    Returns
    -------
    function
        Parameterized IPyWidgets interact function.
    """

    style = {"description_width": "initial"}
    superfeature_ids = [superfeature.id for superfeature in dynophore.superfeatures]

    func = interact(
        plot.static.superfeatures_occurrences,
        dynophore=fixed(dynophore),
        superfeature_names=widgets.SelectMultiple(
            options=["all"] + superfeature_ids,
            value=["all"],
            description="Superfeature name(s):",
            style=style,
        ),
        color_by_feature_type=widgets.Checkbox(value=True, description="Color by feature type"),
        n_equidistant_frames=widgets.IntSlider(
            value=1000,
            min=2,
            max=dynophore.n_frames,
            step=1,
            description="# equidistant frames:",
            style=style,
        ),
    )

    return func


def envpartners_occurrences(dynophore):
    """
    Generate interactive widget to plot the superfeatures' occurrences as barcode.

    Parameters
    ----------
    dynophore : dynophores.Dynophore
        Dynophore.

    Returns
    -------
    function
        Parameterized IPyWidgets interact function.
    """

    style = {"description_width": "initial"}
    superfeature_ids = [superfeature.id for superfeature in dynophore.superfeatures]

    func = interact(
        plot.static.envpartners_occurrences,
        dynophore=fixed(dynophore),
        superfeature_names=widgets.SelectMultiple(
            options=superfeature_ids,
            value=[superfeature_ids[0]],
            description="Superfeature name(s):",
            style=style,
        ),
        n_equidistant_frames=widgets.IntSlider(
            value=1000,
            min=2,
            max=dynophore.n_frames,
            step=1,
            description="# equidistant frames:",
            style=style,
        ),
    )

    return func


def envpartner_distances(dynophore):
    """
    Plot interaction distances for a superfeatures as frame series or histogram.

    Parameters
    ----------
    dynophore : dynophores.Dynophore
        Dynophore.

    Returns
    -------
    function
        Parameterized IPyWidgets interact function.
    """

    style = {"description_width": "initial"}
    superfeature_ids = [superfeature.id for superfeature in dynophore.superfeatures]

    func = interact(
        plot.static.envpartner_distances,
        dynophore=fixed(dynophore),
        superfeature_names=widgets.SelectMultiple(
            options=superfeature_ids,
            value=[superfeature_ids[0]],
            description="Superfeature name(s):",
            style=style,
        ),
        kind=widgets.ToggleButtons(
            options=["line", "hist"],
            description="Plot style",
            button_style="",
            tooltips=["Series", "Histogram"],
        ),
    )

    return func


def envpartners_all_in_one(dynophore):
    """
    Plot interaction data for a superfeature, i.e. occurrences (frame series) and distances
    (frame series and histogram).

    Parameters
    ----------
    dynophore : dynophores.Dynophore
        Dynophore.

    Returns
    -------
    function
        Parameterized IPyWidgets interact function.
    """

    style = {"description_width": "initial"}
    superfeature_ids = [superfeature.id for superfeature in dynophore.superfeatures]

    func = interact(
        plot.static.envpartners_all_in_one,
        dynophore=fixed(dynophore),
        superfeature_name=widgets.Select(
            options=superfeature_ids,
            value=superfeature_ids[0],
            description="Superfeature name(s):",
            style=style,
        ),
        n_equidistant_frames=widgets.IntSlider(
            value=1000,
            min=2,
            max=dynophore.n_frames,
            step=1,
            description="# equidistant frames:",
            style=style,
        ),
    )

    return func
