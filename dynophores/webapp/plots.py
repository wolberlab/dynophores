"""
This module defines plotting functions for the webapp.
"""

import plotly.express as px

from dynophores.plots import _prepare_plot_occurrences


def superfeatures_occurrences_px(dynophore):
    """
    """

    df = dynophore.superfeatures_occurrences
    df.columns.name = "Superfeatures"
    df.index.name = "Frame"
    #df = df.replace(0, None)
    fig = px.scatter(df, facet_row="Superfeatures", facet_row_spacing=0.01)
    # hide and lock down axes
    #fig.update_yaxes(visible=False, fixedrange=True)
    # Remove facet/subplot labels
    fig.update_layout(annotations=[], overwrite=True)
    # strip down the rest of the plot
    fig.update_layout(
        plot_bgcolor="white"
    )
    return fig


def interaction_heatmap(
    dynophore, selector, title="Interaction heatmap", yaxis=""
):
    """
    """

    df = selector(dynophore)

    data_dict = [{"type": "heatmap", "z": df, "name": "Interaction heatmap"}]
    plot_dict = {
        "data": data_dict,
        "layout": {"title": {"text": "Superfeature occurrences"}},
    }
    return plot_dict


def superfeatures_occurrences(
    dynophore, selector, title="Superfeatures occurrences", xaxis="Frames", yaxis="Superfeatures"
):
    """
    Create metadata for plotly boxplot? TODO

    Parameters
    ----------
    dynophore : dynophores.Dynophore
        Dynophore.
    selector : function
        Selects data from dynophore.
    title : str
        Plot title.
    yaxis : str
        Y axis label.
    """

    df = selector(dynophore)
    df = _prepare_plot_occurrences(df)

    data_dict = [
        {"type": "scatter", "x": column.index, "y": column.values, "name": index}
        for (index, column) in df.iteritems()
    ]
    plot_dict = {
        "data": data_dict,
        "layout": {"title": {"text": "Superfeature occurrences"}},
    }
    return plot_dict


def interaction_distances(
    dynophore, selector, title="Interaction distances", xaxis="Frames", yaxis="Distance"
):
    """
    Create metadata for plotly boxplot? TODO

    Parameters
    ----------
    dynophore : dynophores.Dynophore
        Dynophore.
    selector : function
        Selects data from dynophore.
    title : str
        Plot title.
    yaxis : str
        Y axis label.
    """

    df = selector(dynophore)

    data_dict = [
        {"type": "scatter", "x": column.index, "y": column.values, "name": index}
        for (index, column) in df.iteritems()
    ]
    plot_dict = {
        "data": data_dict,
        "layout": {"title": {"text": "Interaction distances"}},
    }
    return plot_dict
