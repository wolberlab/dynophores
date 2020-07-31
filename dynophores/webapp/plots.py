"""
This module defines plotting functions for the webapp.
"""


def boxplot(df, title="Boxplot", yaxis="Error"):
    """
    Create metadata for plotly boxplot? TODO

    Parameters
    ----------
    df : pandas.DataFrame
        TODO
    title : str
        Plot title.
    yaxis : str
        Y axis label.
    """

    data = {"y": df, "type": "box", "name": "blub"}
    data_plot = {
        "data": data,
        "layout": {"title": title, "yaxis": {"title": {"text": yaxis}}},
    }

    return data_plot
