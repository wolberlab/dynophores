"""
Contains static plotting functions.
"""

import math
import itertools

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import ticker
import seaborn as sns

plt.style.use("seaborn")


def superfeatures_vs_envpartners(dynophore, superfeature_ids="all"):
    """
    Plot heatmap of interactions between superfeatures and interaction partners.

    Parameters
    ----------
    dynophore : dynophores.Dynophore
        Dynophore.
    superfeature_ids : str or list of str
        Show all superfeatures (default) or select one or more superfeatures by their superfeature
        identifier.

    Returns
    -------
    fig : matplotlib.figure.Figure
        Plot figure.
    ax : matplotlib.axis.Subplot
        Plot axes.
    """

    superfeature_ids = _format_superfeature_ids(dynophore, superfeature_ids)

    # Prepare data
    data = dynophore.frequency[
        dynophore.frequency.loc["any", :].sort_values(ascending=False).index
    ]
    if superfeature_ids != "all":
        data = data[superfeature_ids]
        data = data[~(data == 0).all(axis=1)]

    fig, ax = plt.subplots(1, 1)
    sns.heatmap(
        data, cmap="Blues", cbar_kws={"label": "Occurrence frequency [%]"}, vmin=0, vmax=100
    )
    ax.set_xlabel("Superfeature IDs")
    ax.set_ylabel("Environmental partner IDs")

    return fig, ax


def superfeatures_occurrences(
    dynophore,
    superfeature_ids="all",
    color_by_feature_type=True,
    frame_range=[0, None],
    frame_step_size=1,
):
    """
    Plot the superfeatures' occurrences as barcode.

    Parameters
    ----------
    dynophore : dynophores.Dynophore
        Dynophore.
    superfeature_ids : str or list of str
        Show all superfeatures (default) or select one or more superfeatures by their superfeature
        identifier.
    color_by_feature_type : bool
        Color barcode by feature type (default) or color all in black.
    frame_range : list of int or list of [int, None]
        Select frame range [start, end]. If end is None, last available frame will be used.
        Default: Select first (0) and last (None) frames.
    frame_step_size : int
        Define frame slicing by step size. Default is 1, i.e. every frame will be selected.
        If e.g. step size is 10, every 10th frame will be selected.

    Returns
    -------
    fig : matplotlib.figure.Figure
        Plot figure.
    ax : matplotlib.axis.Subplot
        Plot axes.
    """

    superfeature_ids = _format_superfeature_ids(dynophore, superfeature_ids)

    # Prepare data
    data = dynophore.superfeatures_occurrences
    if superfeature_ids != "all":
        data = data[superfeature_ids]
    data = _prepare_dataframe_for_plotting(data, frame_range, frame_step_size)

    # If selected data contains no events, do not plot but raise error
    if (data == 0).all().all():
        raise ValueError("Your selection contains no data to plot.")

    # Get frame indices with events
    events_dict = {
        superfeature_id: superfeature[superfeature == 1].index.to_list()
        for superfeature_id, superfeature in data.items()
    }
    events = list(events_dict.values())
    superfeature_ids = list(events_dict.keys())

    # Feature type colors?
    if color_by_feature_type:
        colors = [
            f"#{dynophore.superfeatures[superfeature_id].color}"
            for superfeature_id in data.columns
        ]
    else:
        colors = "black"

    # Plot (plot size depending on number barcodes)
    fig, ax = plt.subplots(figsize=(10, 0.5 + len(events) / 2))
    ax = _occurrences(
        ax=ax,
        events=events,
        colors=colors,
        yticklabels=dynophore._superfeature_ids_frequencies_strings(superfeature_ids),
        xlabel="Frame index",
        xlim=(data.index[0], data.index[-1]),
    )

    return fig, ax


def envpartners_occurrences(
    dynophore,
    superfeature_ids,
    frame_range=[0, None],
    frame_step_size=1,
):
    """
    Plot a superfeature's interaction ocurrences with its interaction partners.

    Parameters
    ----------
    dynophore : dynophores.Dynophore
        Dynophore.
    superfeature_ids : str or list of str
        Superfeature ID(s).
    frame_range : list of int or list of [int, None]
        Select frame range [start, end]. If end is None, last available frame will be used.
        Default: Select first (0) and last (None) frames.
    frame_step_size : int
        Define frame slicing by step size. Default is 1, i.e. every frame will be selected.
        If e.g. step size is 10, every 10th frame will be selected.

    Returns
    -------
    fig : matplotlib.figure.Figure
        Plot figure.
    axes : matplotlib.axis.Subplot or numpy.array of matplotlib.axis.Subplot
        Plot axes.
    """

    superfeature_ids = _format_superfeature_ids(dynophore, superfeature_ids)

    fig, axes = plt.subplots(
        nrows=len(superfeature_ids),
        ncols=1,
        figsize=(10, len(superfeature_ids) * 2),
        sharex=True,
    )

    for i, superfeature_id in enumerate(superfeature_ids):

        if len(superfeature_ids) > 1:
            ax = axes[i]
        else:
            ax = axes

        # Add plot title
        superfeature_title = dynophore._superfeature_ids_frequencies_strings([superfeature_id])
        ax.set_title(superfeature_title[0])

        # Prepare data
        data = dynophore.envpartners_occurrences_by_superfeature(superfeature_id)
        data = _prepare_dataframe_for_plotting(data, frame_range, frame_step_size)
        if data.isna().all().all():
            ax.set_yticks([0])
            ax.set_yticklabels([""])
        else:

            # Get frame indices with events
            events_dict = {
                envpartner_id: envpartner[envpartner == 1].index.to_list()
                for envpartner_id, envpartner in data.items()
            }
            events = list(events_dict.values())
            envpartner_ids = list(events_dict.keys())

            ax = _occurrences(
                ax=ax,
                events=events,
                colors="black",
                yticklabels=dynophore._envpartner_names_frequencies_strings(
                    superfeature_id, envpartner_ids
                ),
                xlabel="",
                xlim=(data.index[0], data.index[-1]),
            )

        if len(superfeature_ids) > 1:
            axes[-1].set_xlabel("Frame index")
        else:
            axes.set_xlabel("Frame index")

    return fig, axes


def _occurrences(ax, events, colors, yticklabels, xlabel, xlim):
    """
    Barcodes for superfeatures and for superfeatures' envpartners are set up in the same way,
    so use this privat helper function to reduce redundant code.

    Parameters
    ----------
    ax : matplotlib.axis.Subplot
        Subplot that shall be used to plot barcode.
    events : list of list
        List of list(s) with frame indices where an event occurs.
    colors : str or int of str
        A single color for all barcodes or colors for each barcode (#colors == #barcodes).
    yticklabels : list of str
        Y axis ticks labels for each barcode (#ylabels == #barcodes).
    xlim : str
        X axis label.

    Returns
    -------
    matplotlib.axis.Subplot
        Subplot with barcode(s).
    """

    ax.eventplot(events, lineoffsets=1, linelength=0.7, linewidths=1, color=colors)
    # Format y axis
    positions = range(0, len(yticklabels))
    # This might be a bug in matplotlib.pyplot.eventplot:
    # If only one barcode is drawn it is centered at 1 (y) instead of 0 (for >1 barcodes).
    # Thus, correct `positions` in case of one barcode from [0] to [1]
    if len(events) == 1:
        positions = [1]  # Instead of [0]
    ax.yaxis.set_major_locator(ticker.FixedLocator(positions))
    ax.yaxis.set_major_formatter(ticker.FixedFormatter(yticklabels))
    # Format x axis
    ax.set_xlabel(xlabel)
    ax.set_xlim(xlim)
    ax.invert_yaxis()

    return ax


def envpartners_distances(
    dynophore, superfeature_ids, kind="line", frame_range=[0, None], frame_step_size=1
):
    """
    Plot interaction distances for a superfeatures as frame series or histogram.

    Parameters
    ----------
    dynophore : dynophores.Dynophore
        Dynophore.
    superfeature_ids : str or list of str
        Show all superfeatures (default) or select one or more superfeatures by their superfeature
        identifier.
    kind : str
        Plot kind, 'line' (distance vs. frames) or 'hist' (distance histogram)
    frame_range : list of int or list of [int, None]
        Select frame range [start, end]. If end is None, last available frame will be used.
        Default: Select first (0) and last (None) frames.
    frame_step_size : int
        Define frame slicing by step size. Default is 1, i.e. every frame will be selected.
        If e.g. step size is 10, every 10th frame will be selected.

    Returns
    -------
    fig : matplotlib.figure.Figure
        Plot figure.
    ax : matplotlib.axis.Subplot
        Plot axes.
    """

    superfeature_ids = _format_superfeature_ids(dynophore, superfeature_ids)

    fig, axes = plt.subplots(
        nrows=len(superfeature_ids),
        ncols=1,
        figsize=(10, len(superfeature_ids) * 4),
        sharex=True,
    )

    for i, superfeature_id in enumerate(superfeature_ids):

        if len(superfeature_ids) > 1:
            ax = axes[i]
        else:
            ax = axes

        # Add plot title
        superfeature_title = dynophore._superfeature_ids_frequencies_strings([superfeature_id])
        ax.set_title(superfeature_title[0])

        # Prepare data
        data = dynophore.envpartners_distances_by_superfeature(superfeature_id)
        data = _prepare_dataframe_for_plotting(data, frame_range, frame_step_size)
        # Add % to environmental partners
        data.columns = dynophore._envpartner_names_frequencies_strings(
            superfeature_id, data.columns.to_list()
        )

        if kind == "line":
            data.plot(kind="line", ax=ax)
            ax.set_xlim((data.index[0], data.index[-1]))
            ax.set_xlabel("Frame index")
            ax.set_ylabel(r"Distance [$\AA$]")
            ax.legend(loc=6, bbox_to_anchor=(1, 0.5), fontsize=12)
        elif kind == "hist":
            value_floor = int(np.floor(data.min().min()))
            value_ceil = int(np.ceil(data.max().max()))
            data.plot(kind="hist", ax=ax, bins=np.arange(value_floor, value_ceil, 0.1), alpha=0.8)
            ax.set_xlim((value_floor, value_ceil))
            ax.set_xlabel(r"Distance [$\AA$]")
            ax.legend(loc=6, bbox_to_anchor=(1, 0.5), fontsize=12)
        else:
            raise KeyError('Plotting kind is unknown. Choose from "line" and "hist".')

    return fig, axes


def envpartners_all_in_one(dynophore, superfeature_id, frame_range=[0, None], frame_step_size=1):
    """
    Plot interaction data for a superfeature, i.e. occurrences (frame series) and distances
    (frame series and histogram).

    Parameters
    ----------
    dynophore : dynophores.Dynophore
        Dynophore.
    superfeature_id : str
        Superfeature ID.
    frame_range : list of int or list of [int, None]
        Select frame range [start, end]. If end is None, last available frame will be used.
        Default: Select first (0) and last (None) frames.
    frame_step_size : int
        Define frame slicing by step size. Default is 1, i.e. every frame will be selected.
        If e.g. step size is 10, every 10th frame will be selected.

    Returns
    -------
    fig : matplotlib.figure.Figure
        Plot figure.
    axes : matplotlib.axis.Subplot
        Plot axes.
    """

    # IPyWidgets' interact function: Cast tuple > str
    if isinstance(superfeature_id, tuple):
        superfeature_id = list(superfeature_id)[0]  # Keep only first
    dynophore._raise_keyerror_if_invalid_superfeature_id(superfeature_id)

    # Get data
    occurrences = dynophore.envpartners_occurrences_by_superfeature(superfeature_id)
    distances = dynophore.envpartners_distances_by_superfeature(superfeature_id)
    occurrences = _prepare_dataframe_for_plotting(occurrences, frame_range, frame_step_size)
    distances = _prepare_dataframe_for_plotting(distances, frame_range, frame_step_size)

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

    """
        # Get frame indices with events
    events_dict = {
        superfeature_id: superfeature[superfeature == 1].index.to_list()
        for superfeature_id, superfeature in data.items()
    }
    events = list(events_dict.values())
    superfeature_ids = list(events_dict.keys())

    # Feature type colors?
    if color_by_feature_type:
        colors = [
            f"#{dynophore.superfeatures[superfeature_id].color}"
            for superfeature_id in data.columns
        ]
    else:
        colors = "black"

    # Plot (plot size depending on number barcodes)
    fig, ax = plt.subplots(figsize=(10, 0.5 + len(events) / 2))
    ax = _occurrences(
        ax=ax,
        events=events,
        colors=colors,
        yticklabels=dynophore._superfeature_ids_frequencies_strings(superfeature_ids),
        xlabel="Frame index",
        xlim=(data.index[0], data.index[-1]),
    )
    """

    # Subplot (0, 1): Empty (will hold legend from subplot (1, 1))
    axes[0][1].axis("off")

    # Subplot (1, 0): Distance time series
    distances.plot(ax=axes[1][0], kind="line", legend=None, linewidth=1)
    axes[1][0].set_xlabel("Frame index", fontsize=16)
    axes[1][0].set_ylabel(r"Distance [$\AA$]", fontsize=16)

    # Subplot (1, 1): Distance histogram
    bins = range(0, math.ceil(distances.max().max()) + 1)
    # Add % to environmental partners
    distances.columns = dynophore._envpartner_names_frequencies_strings(
        superfeature_id, distances.columns.to_list()
    )
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


def _prepare_dataframe_for_plotting(dataframe, frame_range=[0, None], frame_step_size=1):
    """
    Prepare DataFrame for plotting.

    Parameters
    ----------
    dataframe : pandas.DataFrame
        DataFrame of occurrences (or distances) of superfeatures/envpartners (columns) over
        trajectory frames (rows).
    frame_range : list of int or list of [int, None]
        Select frame range [start, end]. If end is None, last available frame will be used.
        Default: Select first (0) and last (None) frames.
    frame_step_size : int
        Define frame slicing by step size. Default is 1, i.e. every frame will be selected.
        If e.g. step size is 10, every 10th frame will be selected.

    Returns
    -------
    pandas.DataFrame
        DataFrame ready for plotting.
    """

    # Sort columns by superfeature/envpartner frequency
    sorted_columns = dataframe.apply(sum).sort_values(ascending=False).index
    dataframe = dataframe[sorted_columns]

    # Slice rows
    dataframe = _slice_dataframe_rows(dataframe, frame_range, frame_step_size)

    return dataframe


def _slice_dataframe_rows(dataframe, ix_range, ix_step_size):
    """
    Slice rows from DataFrame.

    Parameters
    ----------
    data : pandas.DataFrame
        Dataset to be sliced.
    start : int
        Start index.
    end : int
        End index.
    step_size : int
        Step size of indices to be selected.
    """

    ix_start = ix_range[0]
    ix_end = ix_range[1]
    if ix_end is None:
        ix_end = dataframe.shape[0]
    ix_end = ix_end + 1

    selected_indices = itertools.islice(dataframe.index, ix_start, ix_end, ix_step_size)
    dataframe = dataframe.iloc[list(selected_indices), :]

    return dataframe


def _format_superfeature_ids(dynophore, superfeature_ids):
    """
    Format input superfeature IDs (e.g. when they come from IPyWidgets).
    Check if unknown superfeatures were given.

    Parameters
    ----------
    dynophore : dynophores.Dynophore
        Dynophore.
    superfeature_ids : str or list or tuple
        Superfeature IDs.

    Returns
    -------
    str or list of str
        Formated superfeature IDs. String if "all", else list of strings.
    """

    # IPyWidgets' interact function: Cast tuple > list
    if isinstance(superfeature_ids, tuple):
        superfeature_ids = list(superfeature_ids)

    if "all" in superfeature_ids:
        superfeature_ids = "all"
    else:
        if not isinstance(superfeature_ids, list):
            superfeature_ids = [superfeature_ids]
        for superfeature_id in superfeature_ids:
            dynophore._raise_keyerror_if_invalid_superfeature_id(superfeature_id)

    return superfeature_ids
