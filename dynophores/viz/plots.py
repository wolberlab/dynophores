"""
Contains plotting functions for e.g. Jupyter notebooks.
"""

import math

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

from dynophores.definitions import FEATURE_COLORS

########
# Dash #
########


###########
# Jupyter #
###########


def plot_superfeatures_occurrences(
    dynophore, superfeature_names=None, color_by_feature_type=True, max_frames=1000
):
    """
    Plot the superfeatures' occurrences as barcode.

    Parameters
    ----------
    dynophore : dynophores.Dynophore
        Dynophore.
    superfeature_names : None or str or list of str
        Select superfeatures to plot or select all (default).
    color_by_feature_type : bool
        Color barcode by feature type (default) or color all in black.
    max_frames : int
        Number of frames to display in barcode plot. If input data contains more than `max_frames`,
        `max_frames` equidistant frames will be selected.
    """

    occurrences = dynophore.superfeatures_occurrences

    if superfeature_names is not None:
        if isinstance(superfeature_names, str):
            superfeature_names = [superfeature_names]

        # Get all superfeature names that are in data
        superfeature_names_curated = [i for i in superfeature_names if i in occurrences.columns]
        superfeature_names_omitted = list(
            set(superfeature_names) - set(superfeature_names_curated)
        )
        if len(superfeature_names_omitted) > 0:
            print(f"Superfeature names {superfeature_names_omitted} omitted because unknown.")

        # Select subset
        occurrences = occurrences[superfeature_names_curated]

    # Prepare data
    occurrences = _prepare_plot_occurrences(occurrences, max_frames)

    # Feature type colors?
    if color_by_feature_type:
        feature_types = [i.split("[")[0] for i in occurrences.columns]
        colors = [
            FEATURE_COLORS[i] if i in FEATURE_COLORS.keys() else "black" for i in feature_types
        ]
    else:
        colors = None

    # Plot (plot size depending on number barcodes)
    fig, ax = plt.subplots(figsize=(10, occurrences.shape[1] / 2))
    occurrences.plot(marker=".", markersize=5, linestyle="", legend=None, ax=ax, color=colors)
    # Set y tick labels
    ax.set_yticks(range(0, occurrences.shape[1] + 2))
    ax.set_yticklabels([""] + occurrences.columns.to_list() + [""])
    ax.invert_yaxis()
    # Set x axis limits and label
    ax.set_xlabel("frame index")
    ax.set_xlim((occurrences.index[0], occurrences.index[-1]))


def plot_envpartners_occurrences(dynophore, superfeature_name, max_frames=1000):
    """
    Plot a superfeature's interaction ocurrences with its interaction partners.

    Parameters
    ----------
    dynophore : dynophores.Dynophore
        Dynophore.
    superfeature_name : str
        Superfeature name
    max_frames : int
        Number of frames to display in barcode plot. If input data contains more than `max_frames`,
        `max_frames` equidistant frames will be selected.
    """

    dynophore.is_superfeature(superfeature_name)

    # Prepare data
    occurrences = _prepare_plot_occurrences(
        dynophore.envpartners_occurrences[superfeature_name], max_frames
    )

    # Plot (plot size depending on number barcodes)
    fig, ax = plt.subplots(figsize=(10, occurrences.shape[1] / 2))
    occurrences.plot(marker=".", markersize=5, linestyle="", legend=None, ax=ax)
    # Set y tick labels
    ax.set_yticks(range(0, occurrences.shape[1] + 2))
    ax.set_yticklabels([""] + occurrences.columns.to_list() + [""])
    ax.invert_yaxis()
    # Set x axis limits and label
    ax.set_xlabel("Frame")
    ax.set_xlim((occurrences.index[0], occurrences.index[-1]))


def plot_envpartners(dynophore, superfeature_name, max_frames=1000):
    """
    Plot interaction data for a superfeature, i.e. occurrences (frame series) and distances
    (frame series and histogram).

    Parameters
    ----------
    dynophore : dynophores.Dynophore
        Dynophore.
    superfeature_name : str
        Superfeature name
    max_frames : int
        Number of frames to display in barcode plot. If input data contains more than `max_frames`,
        `max_frames` equidistant frames will be selected.
    """

    dynophore.is_superfeature(superfeature_name)
    occurrences = _prepare_plot_envparters_occurrences(dynophore, superfeature_name, max_frames)
    distances = _prepare_plot_envpartners_distances(dynophore, superfeature_name, max_frames)

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
    axes[0][0].set_xlabel("frame index")
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
        alpha=0.7,
        density=True,
        xlim=(0, 1),
    )
    axes[1][1].set_xlabel("Frequency", fontsize=16)
    axes[1][1].legend(loc=6, bbox_to_anchor=(0, 1.5), fontsize=12)


def plot_superfeatures_vs_envpartners(dynophore):
    """
    Plot heatmap of interactions between superfeatures and interaction partners.

    Parameters
    ----------
    dynophore : dynophores.Dynophore
        Dynophore.
    """

    # Sort superfeatures by overall frequency
    data = dynophore.frequency[
        dynophore.frequency.loc["any", :].sort_values(ascending=False).index
    ]
    data.fillna(0, inplace=True)
    sns.heatmap(data, cmap="Blues")


def plot_envpartner_distances(dynophore, superfeature_name, kind):
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
    """

    data = dynophore.envpartners_distances[superfeature_name]

    if kind == "line":
        fig, ax = plt.subplots(figsize=(10, 5))
        data.plot(kind="line", ax=ax)
        ax.set_xlim((0, data.shape[0]))
        ax.set_xlabel("Frame index")
        ax.set_ylabel(r"Distance [$\AA$]")
    elif kind == "hist":
        ax = data.plot(kind="hist")
        ax.set_xlabel(r"Distance [$\AA$]")
    else:
        raise ValueError('Plotting kind is unknown. Choose from "line" and "hist".')


def _prepare_plot_occurrences(occurrences, max_frames=1000):
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
    max_frames : int
        Number of frames to display in barcode plot. If input data contains more than `max_frames`,
        `max_frames` equidistant frames will be selected.

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
        selected_indices = [i for i in range(0, 1000, math.floor(1002 / max_frames))]
        occurrences = occurrences.iloc[selected_indices, :]

    # Transform 1 in binary values to rank in plot
    occurrences_plot = {}
    for i, (name, data) in enumerate(occurrences.items()):
        data = data.replace([0, 1], [None, i + 1])
        occurrences_plot[name] = data
    occurrences_plot = pd.DataFrame(occurrences_plot)

    return occurrences_plot


def _prepare_plot_envparters_occurrences(dynophore, superfeature_name, max_frames=1000):
    """
    Prepare data for plotting superfeature interaction occurrences.

    Parameters
    ----------
    dynophore : dynophores.Dynophore
        Dynophore.
    superfeature_name : str
        Superfeature name
    max_frames : int
        Number of frames to display in barcode plot. If input data contains more than `max_frames`,
        `max_frames` equidistant frames will be selected.

    Returns
    -------
    pandas.DataFrame
        Occurrences, ready for plotting.
    """

    occurrences = dynophore.envpartners_occurrences[superfeature_name]
    occurrences = _prepare_plot_occurrences(occurrences, max_frames)

    return occurrences


def _prepare_plot_envpartners_distances(dynophore, superfeature_name, max_frames=1000):
    """
    Prepare data for plotting interaction distances for a superfeature:
    - Sort by interaction frequency

    Parameters
    ----------
    dynophore : dynophores.Dynophore
        Dynophore.
    superfeature_name : str
        Superfeature name
    max_frames : int
        Number of frames to display in barcode plot. If input data contains more than `max_frames`,
        `max_frames` equidistant frames will be selected.

    Returns
    -------
    pandas.DataFrame
        Interaction distances for a superfeature, ready for plotting.
    """

    # Get data
    distances = dynophore.envpartners_distances[superfeature_name]

    # Sort data by ratio
    ratio = (
        dynophore.frequency.loc[:, superfeature_name]
        .dropna()
        .drop("any")
        .sort_values(ascending=False)
    )
    distances = distances[ratio.index]

    return distances
