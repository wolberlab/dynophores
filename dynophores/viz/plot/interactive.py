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

import math

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
    """

    style = {"description_width": "initial"}
    superfeature_ids = (
        dynophore.frequency.loc["any", :].sort_values(ascending=False).index.to_list()
    )
    superfeature_ids = ["all"] + superfeature_ids

    interact(
        plot.static.superfeatures_vs_envpartners,
        dynophore=fixed(dynophore),
        superfeature_ids=widgets.SelectMultiple(
            options=superfeature_ids,
            value=["all"],
            description="Superfeature ID(s):",
            style=style,
        ),
        annotate_heatmap=widgets.Checkbox(value=False, description="Annotate heatmap cells"),
    )


def superfeatures_occurrences(dynophore):
    """
    Generate interactive widget to plot the superfeatures' occurrences as barcode.

    Parameters
    ----------
    dynophore : dynophores.Dynophore
        Dynophore.
    """

    style = {"description_width": "initial"}
    superfeature_ids = (
        dynophore.frequency.loc["any", :].sort_values(ascending=False).index.to_list()
    )

    interact(
        plot.static.superfeatures_occurrences,
        dynophore=fixed(dynophore),
        superfeature_ids=widgets.SelectMultiple(
            options=["all"] + superfeature_ids,
            value=["all"],
            description="Superfeature ID(s):",
            style=style,
        ),
        color_by_feature_type=widgets.Checkbox(value=True, description="Color by feature type"),
        frame_range=widgets.IntRangeSlider(
            value=[0, dynophore.n_frames],
            min=0,
            max=dynophore.n_frames,
            step=1,
            description="Frames range:",
            style=style,
        ),
        frame_step_size=widgets.BoundedIntText(
            # Default value results in displayed 1000 frames
            value=math.floor(dynophore.n_frames / 1000),
            min=1,
            max=dynophore.n_frames - 1,
            step=1,
            description="Frames step size:",
            style=style,
        ),
    )


def envpartners_occurrences(dynophore):
    """
    Generate interactive widget to plot the superfeatures' occurrences as barcode.

    Parameters
    ----------
    dynophore : dynophores.Dynophore
        Dynophore.
    """

    style = {"description_width": "initial"}
    superfeature_ids = (
        dynophore.frequency.loc["any", :].sort_values(ascending=False).index.to_list()
    )

    interact(
        plot.static.envpartners_occurrences,
        dynophore=fixed(dynophore),
        superfeature_ids=widgets.SelectMultiple(
            options=superfeature_ids,
            value=[superfeature_ids[0]],
            description="Superfeature ID(s):",
            style=style,
        ),
        frame_range=widgets.IntRangeSlider(
            value=[0, dynophore.n_frames],
            min=0,
            max=dynophore.n_frames,
            step=1,
            description="Frames range:",
            style=style,
        ),
        frame_step_size=widgets.BoundedIntText(
            # Default value results in displayed 1000 frames
            value=math.floor(dynophore.n_frames / 1000),
            min=1,
            max=dynophore.n_frames - 1,
            step=1,
            description="Frames step size:",
            style=style,
        ),
        occurrence_min=widgets.BoundedFloatText(
            value=0,
            min=0,
            max=100,
            step=1,
            description="Occurrence minimum [%]:",
            style=style,
        ),
    )


def envpartners_distances(dynophore):
    """
    Plot interaction distances for a superfeatures as frame series or histogram.

    Parameters
    ----------
    dynophore : dynophores.Dynophore
        Dynophore.
    """

    style = {"description_width": "initial"}
    superfeature_ids = (
        dynophore.frequency.loc["any", :].sort_values(ascending=False).index.to_list()
    )

    interact(
        plot.static.envpartners_distances,
        dynophore=fixed(dynophore),
        superfeature_ids=widgets.SelectMultiple(
            options=superfeature_ids,
            value=[superfeature_ids[0]],
            description="Superfeature ID(s):",
            style=style,
        ),
        kind=widgets.ToggleButtons(
            options=["line", "hist"],
            description="Plot style",
            button_style="",
            tooltips=["Series", "Histogram"],
        ),
        frame_range=widgets.IntRangeSlider(
            value=[0, dynophore.n_frames],
            min=0,
            max=dynophore.n_frames,
            step=1,
            description="Frames range:",
            style=style,
        ),
        frame_step_size=widgets.BoundedIntText(
            # Default value results in displayed 1000 frames
            value=math.floor(dynophore.n_frames / 1000),
            min=1,
            max=dynophore.n_frames - 1,
            step=1,
            description="Frames step size:",
            style=style,
        ),
        occurrence_min=widgets.BoundedFloatText(
            value=0,
            min=0,
            max=100,
            step=1,
            description="Occurrence minimum [%]:",
            style=style,
        ),
    )


def envpartners_all_in_one(dynophore):
    """
    Plot interaction data for a superfeature, i.e. occurrences (frame series) and distances
    (frame series and histogram).

    Parameters
    ----------
    dynophore : dynophores.Dynophore
        Dynophore.
    """

    style = {"description_width": "initial"}
    superfeature_ids = (
        dynophore.frequency.loc["any", :].sort_values(ascending=False).index.to_list()
    )

    interact(
        plot.static.envpartners_all_in_one,
        dynophore=fixed(dynophore),
        superfeature_id=widgets.Select(
            options=superfeature_ids,
            value=superfeature_ids[0],
            description="Superfeature ID(s):",
            style=style,
        ),
        frame_range=widgets.IntRangeSlider(
            value=[0, dynophore.n_frames],
            min=0,
            max=dynophore.n_frames,
            step=1,
            description="Frames range:",
            style=style,
        ),
        frame_step_size=widgets.BoundedIntText(
            # Default value results in displayed 1000 frames
            value=math.floor(dynophore.n_frames / 1000),
            min=1,
            max=dynophore.n_frames - 1,
            step=1,
            description="Frames step size:",
            style=style,
        ),
        occurrence_min=widgets.BoundedFloatText(
            value=0,
            min=0,
            max=100,
            step=1,
            description="Occurrence minimum [%]:",
            style=style,
        ),
    )
