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


def superfeatures_vs_envpartners(dynophore, superfeature_ids="all", annotate_heatmap=False):
    """
    Plot heatmap of interactions between superfeatures and interaction partners.

    Parameters
    ----------
    dynophore : dynophores.Dynophore
        Dynophore.
    superfeature_ids : str or list of str
        Show all superfeatures (default) or select one or more superfeatures by their superfeature
        identifier.
    annotate_heatmap : bool
        Annotate heatmap cells (default: False).

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

    fig, ax = plt.subplots(1, 1, tight_layout=True)
    sns.heatmap(
        data,
        annot=annotate_heatmap,
        cmap="Blues",
        cbar_kws={"label": "Occurrence frequency [%]"},
        vmin=0,
        vmax=100,
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
    data = _prepare_superfeature_plotting_data(
        dynophore, superfeature_ids, frame_range, frame_step_size
    )

    # If selected data contains no events, do not plot but raise error
    if (data == 0).all().all() or data.shape[1] == 0:
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
    fig, ax = plt.subplots(figsize=(8, 0.5 + len(events) / 2), tight_layout=True)
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
    dynophore, superfeature_ids, frame_range=[0, None], frame_step_size=1, occurrence_min=0
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
    occurrence_min : int or float
        Remove all envpartners below the occurrence cutoff (default: 0).

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
        figsize=(8, len(superfeature_ids) * 2),
        sharex=True,
        tight_layout=True,
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
        data, _ = _prepare_envpartner_plotting_data(
            dynophore, superfeature_id, frame_range, frame_step_size, occurrence_min
        )
        if (data == 0).all().all() or data.shape[1] == 0:
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


def envpartners_distances(
    dynophore,
    superfeature_ids,
    kind="line",
    frame_range=[0, None],
    frame_step_size=1,
    occurrence_min=0,
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
    occurrence_min : int or float
        Remove all envpartners below the occurrence cutoff (default: 0).

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
        figsize=(8, len(superfeature_ids) * 4),
        sharex=True,
        tight_layout=True,
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
        occurrences, distances = _prepare_envpartner_plotting_data(
            dynophore, superfeature_id, frame_range, frame_step_size, occurrence_min
        )
        # Use all distances
        data = distances
        # Use interaction-relevant distances only:
        # Where occurrence is 0, distances will be 0; set those distances to None
        data_interaction_frames = (distances * occurrences).applymap(
            lambda x: x if x > 0 else None
        )

        if data.shape[1] == 0:
            ax.set_yticks([0])
            ax.set_yticklabels([""])
        else:

            # Add % to environmental partners
            data.columns = dynophore._envpartner_names_frequencies_strings(
                superfeature_id, data.columns.to_list()
            )

            if kind == "line":

                # Plot all distances as line plot
                data.plot(kind="line", ax=ax, linewidth=0.5)

                # Plot all interaction-relevant distances as dot plot
                data_interaction_frames.plot(ax=ax, style=".", markersize=5)

                ax.set_xlim((data.index[0], data.index[-1]))
                ax.set_xlabel("Frame index")
                ax.set_ylabel(r"Distance [$\AA$]")
                # Currently we have legends from both plots; we only need it once!
                legend_labels = data.columns.to_list()
                ax.legend(labels=legend_labels, loc=6, bbox_to_anchor=(1, 0.5))

            elif kind == "hist":

                # TODO or use min/max in full dataset instead of per superfeature?
                value_floor = int(np.floor(data_interaction_frames.min().min()))
                value_ceil = int(np.ceil(data_interaction_frames.max().max()))

                data_interaction_frames.plot(
                    kind="hist", ax=ax, bins=np.arange(value_floor, value_ceil, 0.1), alpha=0.8
                )

                ax.set_xlim((value_floor, value_ceil))
                ax.set_xlabel(r"Distance [$\AA$]")
                ax.legend(loc=6, bbox_to_anchor=(1, 0.5))

            else:
                raise KeyError("Plotting kind is unknown. Choose from 'line' and 'hist'.")

    return fig, axes


def envpartners_all_in_one(
    dynophore, superfeature_id, frame_range=[0, None], frame_step_size=1, occurrence_min=0
):
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
    occurrence_min : int or float
        Remove all envpartners below the occurrence cutoff (default: 0).

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
    occurrences, distances = _prepare_envpartner_plotting_data(
        dynophore, superfeature_id, frame_range, frame_step_size, occurrence_min
    )
    # If selected data contains no events, do not plot but raise error
    if (occurrences.shape[1] == 0) or (distances.shape[1] == 0):
        raise ValueError("Your selection contains no data to plot.")

    # Set up plot
    fig, axes = plt.subplots(
        nrows=2,
        ncols=2,
        figsize=(8, 5),
        sharey="row",
        sharex="col",
        gridspec_kw={"width_ratios": [3, 1], "wspace": 0.05, "hspace": 0.15},
        constrained_layout=True,
    )
    fig.suptitle(superfeature_id)

    # Subplot (0, 0): Interaction occurrences (barplot)
    events_dict = {
        envpartner_id: envpartner[envpartner == 1].index.to_list()
        for envpartner_id, envpartner in occurrences.items()
    }
    events = list(events_dict.values())
    # This is an unpretty fix: eventplot does not circle over default colors, so we have to get
    # default colors and add a few circles to the list
    colors = plt.rcParams["axes.prop_cycle"].by_key()["color"] * 5
    envpartner_ids = list(events_dict.keys())
    axes[0][0] = _occurrences(
        ax=axes[0][0],
        events=events,
        colors=colors[: len(envpartner_ids)],
        yticklabels=dynophore._envpartner_names_frequencies_strings(
            superfeature_id, envpartner_ids
        ),
        xlabel="",
        xlim=(occurrences.index[0], occurrences.index[-1]),
    )

    # Subplot (0, 1): Empty (will hold legend from subplot (1, 1))
    axes[0][1].axis("off")

    # Subplot (1, 0): Distance time series
    distances.plot(ax=axes[1][0], kind="line", legend=None, linewidth=1)
    axes[1][0].set_xlabel("Frame index")
    axes[1][0].set_ylabel(r"Distance [$\AA$]")

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
        legend=None,
    )
    axes[1][1].set_xlabel("Frequency")

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


def _prepare_envpartner_plotting_data(
    dynophore,
    superfeature_id,
    frame_range=[0, None],
    frame_step_size=1,
    occurrence_min=0,
):

    """
    Prepare envpartner data for plotting.

    Parameters
    ----------
    dynophore : dynophores.Dynophore
        Dynophore.
    superfeature_ids : list of str or str
        Superfeature ID.
    frame_range : list of int or list of [int, None]
        Select frame range [start, end]. If end is None, last available frame will be used.
        Default: Select first (0) and last (None) frames.
    frame_step_size : int
        Define frame slicing by step size. Default is 1, i.e. every frame will be selected.
        If e.g. step size is 10, every 10th frame will be selected.
    occurrence_min : int or float
        Remove all envpartners below the occurrence cutoff (default: 0).

    Returns
    -------
    pandas.DataFrame
        DataFrame ready for plotting.
    """

    superfeature = dynophore.superfeatures[superfeature_id]

    # Prepare occurrences
    occurrences = superfeature.envpartners_occurrences
    # Sort columns by envpartner frequency
    sorted_columns = occurrences.sum().sort_values(ascending=False).index
    occurrences = occurrences[sorted_columns]
    # Select columns by occurrence_min
    frequency = superfeature.frequency.drop("any")
    selected_columns = frequency[frequency >= occurrence_min].index
    occurrences = occurrences[selected_columns]
    # Slice rows
    occurrences = _slice_dataframe_rows(occurrences, frame_range, frame_step_size)

    # Prepare distances
    distances = dynophore.envpartners_distances_by_superfeature(superfeature_id)
    # Select distances by occurrences columns/indices
    distances = distances.loc[occurrences.index, occurrences.columns]

    return occurrences, distances


def _prepare_superfeature_plotting_data(
    dynophore, superfeature_ids, frame_range=[0, None], frame_step_size=1
):
    """
    Prepare superfeature data for plotting.

    Parameters
    ----------
    dynophore : dynophores.Dynophore
        Dynophore.
    superfeature_ids : list of str or str
        List of superfeature IDs or "all".
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

    data = dynophore.superfeatures_occurrences
    if superfeature_ids != "all":
        data = data[superfeature_ids]

    # Sort columns by superfeature frequency
    sorted_columns = data.sum().sort_values(ascending=False).index
    data = data[sorted_columns]

    # Select rows
    data = _slice_dataframe_rows(data, frame_range, frame_step_size)

    return data


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
