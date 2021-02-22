"""
Contains static and interactive plotting functions for Jupyter notebooks.

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

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from ipywidgets import interact, fixed
import ipywidgets as widgets

from dynophores.definitions import FEATURE_COLORS

plt.style.use("seaborn")


def superfeatures_occurrences(
    dynophore, superfeature_names="all", color_by_feature_type=True, n_equidistant_frames=1000
):
    """
    Plot the superfeatures' occurrences as barcode.

    Parameters
    ----------
    dynophore : dynophores.Dynophore
        Dynophore.
    superfeature_names : str or list of str
        Show all superfeatures (default) or select one or more superfeatures by their superfeature
        identifier.
    color_by_feature_type : bool
        Color barcode by feature type (default) or color all in black.
    n_equidistant_frames : int
        Number of frames to display in barcode plot. If input data contains more than
        `n_equidistant_frames`, `n_equidistant_frames` equidistant frames will be selected.

    Returns
    -------
    fig : matplotlib.figure.Figure
        Plot figure.
    ax : matplotlib.axis.Subplot
        Plot axes.
    """

    occurrences = dynophore.superfeatures_occurrences

    # Correct input from ipywidgets' interact function
    if isinstance(superfeature_names, tuple):
        superfeature_names = list(superfeature_names)
    if "all" in superfeature_names:
        superfeature_names = "all"

    # (Optionally) select superfeature subset
    if superfeature_names != "all":
        if not isinstance(superfeature_names, list):
            superfeature_names = [superfeature_names]
        for superfeature_name in superfeature_names:
            dynophore.raise_keyerror_if_invalid_superfeature_name(superfeature_name)
        occurrences = occurrences[superfeature_names]

    # Prepare data
    occurrences = _prepare_plot_occurrences(occurrences, n_equidistant_frames)

    # Feature type colors?
    if color_by_feature_type:
        feature_types = [i.split("[")[0] for i in occurrences.columns]
        colors = [
            FEATURE_COLORS[i] if i in FEATURE_COLORS.keys() else "black" for i in feature_types
        ]
    else:
        colors = "black"

    # Plot (plot size depending on number barcodes)
    fig, ax = plt.subplots(figsize=(10, occurrences.shape[1] / 2))
    occurrences.plot(marker=".", markersize=5, linestyle="", legend=None, ax=ax, color=colors)
    # Set y tick labels
    ax.set_yticks(range(0, occurrences.shape[1] + 2))
    ax.set_yticklabels([""] + occurrences.columns.to_list() + [""])
    ax.invert_yaxis()
    # Set x axis limits and label
    ax.set_xlabel("Frame index")
    ax.set_xlim((occurrences.index[0], occurrences.index[-1]))

    return fig, ax


def superfeatures_occurrences_interactive(dynophore):
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
        superfeatures_occurrences,
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


def envpartners_occurrences(dynophore, superfeature_names, n_equidistant_frames=1000):
    """
    Plot a superfeature's interaction ocurrences with its interaction partners.

    Parameters
    ----------
    dynophore : dynophores.Dynophore
        Dynophore.
    superfeature_names : str
        Superfeature name
    n_equidistant_frames : int
        Number of frames to display in barcode plot. If input data contains more than
        `n_equidistant_frames`, `n_equidistant_frames` equidistant frames will be selected.

    Returns
    -------
    fig : matplotlib.figure.Figure
        Plot figure.
    ax : matplotlib.axis.Subplot
        Plot axes.
    """

    # Correct input from ipywidgets' interact function
    if isinstance(superfeature_names, tuple):
        superfeature_names = list(superfeature_names)
    if isinstance(superfeature_names, str):
        superfeature_names = [superfeature_names]
    for superfeature_name in superfeature_names:
        dynophore.raise_keyerror_if_invalid_superfeature_name(superfeature_name)

    # Plot (plot size depending on number barcodes)
    fig, axes = plt.subplots(nrows=len(superfeature_names), ncols=1, sharex=True)
    # figsize=(10, occurrences.shape[1] / 2)

    for i, superfeature_name in enumerate(superfeature_names):

        if len(list(axes)) > 1:
            ax = axes[i]
        else:
            ax = axes

        # Prepare data
        occurrences = _prepare_plot_occurrences(
            dynophore.envpartners_occurrences[superfeature_name], n_equidistant_frames
        )
        occurrences.plot(marker=".", markersize=5, linestyle="", legend=None, ax=ax, color="black")
        # Set y tick labels
        ax.set_yticks(range(0, occurrences.shape[1] + 2))
        ax.set_yticklabels([""] + occurrences.columns.to_list() + [""])
        ax.invert_yaxis()
        # Set x axis limits and label
        ax.set_xlabel("Frame")
        ax.set_xlim((occurrences.index[0], occurrences.index[-1]))

    return fig, axes


def envpartners_occurrences_interactive(dynophore):
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
        envpartners_occurrences,
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


def envpartners(dynophore, superfeature_name, n_equidistant_frames=1000):
    """
    Plot interaction data for a superfeature, i.e. occurrences (frame series) and distances
    (frame series and histogram).

    Parameters
    ----------
    dynophore : dynophores.Dynophore
        Dynophore.
    superfeature_name : str
        Superfeature name
    n_equidistant_frames : int
        Number of frames to display in barcode plot. If input data contains more than
        `n_equidistant_frames`, `n_equidistant_frames` equidistant frames will be selected.

    Returns
    -------
    fig : matplotlib.figure.Figure
        Plot figure.
    axes : matplotlib.axis.Subplot
        Plot axes.
    """

    dynophore.raise_keyerror_if_invalid_superfeature_name(superfeature_name)
    occurrences = _prepare_plot_envparters_occurrences(
        dynophore, superfeature_name, n_equidistant_frames
    )
    distances = _prepare_plot_envpartners_distances(
        dynophore, superfeature_name, n_equidistant_frames
    )

    # Set up plot
    fig, axes = plt.subplots(
        nrows=2,
        ncols=2,
        figsize=(15, 7),
        sharey="row",
        sharex="col",
        gridspec_kw={"width_ratios": [3, 1], "wspace": 0.05, "hspace": 0.05},
    )

    # Subplot (0, 0): Interaction occurrences (barplot)
    occurrences.plot(ax=axes[0][0], kind="line", legend=None, marker=".", linestyle="")
    # Set y tick labels (tick per envpartner but do not show label)
    axes[0][0].set_yticks(range(0, occurrences.shape[1] + 2))
    axes[0][0].set_yticklabels("")
    # Set x axis limits and label
    axes[0][0].set_xlabel("Frame index")
    axes[0][0].set_xlim((occurrences.index[0], occurrences.index[-1]))
    # Show interactions from top to bottom from most common to rarest interactions
    axes[0][0].invert_yaxis()

    # Subplot (0, 1): Empty (will hold legend from subplot (1, 1))
    axes[0][1].axis("off")

    # Subplot (1, 0): Distance time series
    distances.plot(ax=axes[1][0], kind="line", legend=None, linewidth=1)
    axes[1][0].set_xlabel("Frame index", fontsize=16)
    axes[1][0].set_ylabel(r"Distance [$\AA$]", fontsize=16)

    # Subplot (1, 1): Distance histogram
    bins = range(0, math.ceil(distances.max().max()) + 1)
    distances.plot(
        ax=axes[1][1],
        kind="hist",
        bins=bins,
        orientation="horizontal",
        alpha=0.8,
        density=True,
        xlim=(0, 1),
    )
    axes[1][1].set_xlabel("Frequency", fontsize=16)
    axes[1][1].legend(loc=6, bbox_to_anchor=(0, 1.5), fontsize=12)

    return fig, axes


def superfeatures_vs_envpartners(dynophore):
    """
    Plot heatmap of interactions between superfeatures and interaction partners.

    Parameters
    ----------
    dynophore : dynophores.Dynophore
        Dynophore.

    Returns
    -------
    fig : matplotlib.figure.Figure
        Plot figure.
    ax : matplotlib.axis.Subplot
        Plot axes.
    """

    # Sort superfeatures by overall frequency
    data = dynophore.frequency[
        dynophore.frequency.loc["any", :].sort_values(ascending=False).index
    ]
    data.fillna(0, inplace=True)
    fig, ax = plt.subplots(1, 1)
    sns.heatmap(data, cmap="Blues")

    return fig, ax


def envpartner_distances(dynophore, superfeature_name, kind="line"):
    """
    Plot interaction distances for a superfeatures as frame series or histogram.

    Parameters
    ----------
    dynophore : dynophores.Dynophore
        Dynophore.
    superfeature_name : str
        Superfeature name.
    kind : str
        Plot kind, 'line' (distance vs. frames) or 'hist' (distance histogram)

    Returns
    -------
    fig : matplotlib.figure.Figure
        Plot figure.
    ax : matplotlib.axis.Subplot
        Plot axes.
    """

    dynophore.raise_keyerror_if_invalid_superfeature_name(superfeature_name)

    data = dynophore.envpartners_distances[superfeature_name]

    fig, ax = plt.subplots(figsize=(10, 5))

    if kind == "line":
        data.plot(kind="line", ax=ax)
        ax.set_xlim((0, data.shape[0]))
        ax.set_xlabel("Frame index")
        ax.set_ylabel(r"Distance [$\AA$]")
    elif kind == "hist":
        value_floor = int(np.floor(data.min().min()))
        value_ceil = int(np.ceil(data.max().max()))
        data.plot(kind="hist", ax=ax, bins=np.arange(value_floor, value_ceil, 0.1), alpha=0.8)
        ax.set_xlim((value_floor, value_ceil))
        ax.set_xlabel(r"Distance [$\AA$]")
    else:
        raise KeyError('Plotting kind is unknown. Choose from "line" and "hist".')

    return fig, ax


def _prepare_plot_occurrences(occurrences, n_equidistant_frame=1000):
    """
    Prepare data for plotting occurrences (superfeatures or superfeature interactions):
    - Sort by interaction frequency
    - Get subset of data if more than 1000 frames
    - Transform 1 in binary values to rank in plot

    Parameters
    ----------
    occurrences : pandas.DataFrame
        Occurrences (0 or 1) per frame (rows) for one or more superfeatures or superfeature
        interactions (columns).
    n_equidistant_frames : int
        Number of frames to display in barcode plot.
        If input data contains more than `n_equidistant_frames`, `n_equidistant_frames` equidistant
        frames will be selected.

    Returns
    -------
    pandas.DataFrame
        Occurrences, ready for plotting.
    """

    # Sort data by ratio
    ratio = round(occurrences.apply(sum) / occurrences.shape[0] * 100, 2).sort_values(
        ascending=False
    )
    occurrences = occurrences[ratio.index]

    # Get subset of data if more than 1000 frames
    if occurrences.shape[0] > 1000:
        selected_indices = [i for i in range(0, 1000, math.floor(1002 / n_equidistant_frame))]
        occurrences = occurrences.iloc[selected_indices, :]

    # Transform 1 in binary values to rank in plot
    occurrences_plot = {}
    for i, (name, data) in enumerate(occurrences.items()):
        data = data.replace([0, 1], [None, i + 1])
        occurrences_plot[name] = data
    occurrences_plot = pd.DataFrame(occurrences_plot)

    return occurrences_plot


def _prepare_plot_envparters_occurrences(dynophore, superfeature_name, n_equidistant_frames=1000):
    """
    Prepare data for plotting superfeature interaction occurrences.

    Parameters
    ----------
    dynophore : dynophores.Dynophore
        Dynophore.
    superfeature_name : str
        Superfeature name
    n_equidistant_frames : int
        Number of frames to display in barcode plot. If input data contains more than
        `n_equidistant_frames`, `n_equidistant_frames` equidistant frames will be selected.

    Returns
    -------
    pandas.DataFrame
        Occurrences, ready for plotting.
    """

    occurrences = dynophore.envpartners_occurrences[superfeature_name]
    occurrences = _prepare_plot_occurrences(occurrences, n_equidistant_frames)

    return occurrences


def _prepare_plot_envpartners_distances(dynophore, superfeature_name, n_equidistant_frame=1000):
    """
    Prepare data for plotting interaction distances for a superfeature:
    - Sort by interaction frequency

    Parameters
    ----------
    dynophore : dynophores.Dynophore
        Dynophore.
    superfeature_name : str
        Superfeature name
    n_equidistant_frames : int
        Number of frames to display in barcode plot. If input data contains more than
        `n_equidistant_frames`, `n_equidistant_frames` equidistant frames will be selected.

    Returns
    -------
    pandas.DataFrame
        Interaction distances for a superfeature, ready for plotting.
    """

    # Get data
    distances = dynophore.envpartners_distances[superfeature_name]

    # Sort data by frequency
    frequency = dynophore.frequency[superfeature_name]
    frequency = frequency[frequency > 0].drop("any").sort_values(ascending=False)
    distances = distances[frequency.index]

    return distances
