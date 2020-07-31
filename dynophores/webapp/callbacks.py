"""
This module defines the webapp callbacks.
"""

import base64
import io
import logging
import traceback

import dash
from dash.dependencies import Input, Output, State
import dash_html_components as html
import dash_table
import pandas as pd

from .layout import PLOTS, CARDS
from .content import load_dynophore


from .plots import interaction_distances

logger = logging.getLogger(__name__)


def parse_contents(contents, filename):
    """
    Parse contents from file upload.

    Parameters
    ----------
    contents : str
        File content
    filename : str
        File name.

    Returns
    -------
    dash_html_components.Div.Div
        Output HTML code.
    """

    content_type, content_string = contents.split(',')

    decoded = base64.b64decode(content_string)

    try:
        if 'csv' in filename:
            # Assume that the user uploaded a CSV file
            df = pd.read_csv(
                io.StringIO(decoded.decode('utf-8')))
        elif 'xls' in filename:
            # Assume that the user uploaded an excel file
            df = pd.read_excel(io.BytesIO(decoded))
        else:
            return html.Div([
                f"Selected file: {filename}: This file's format is not supported for upload."
            ])
    except Exception as e:
        return html.Div([
            f"Selected file: {filename}: There was an error processing this file: {e}."
        ])

    html_result = html.Div([
        html.Div(f"File uploaded: {filename}"),
        dash_table.DataTable(
            data=df.to_dict('records'),
            columns=[{'name': i, 'id': i} for i in df.columns]
        ),
    ])

    return html_result


def register_callbacks(app):
    """
    Register callbacks.

    Parameters
    ----------
    app : TODO
        TODO
    """

    figure_outputs = [
        [
            Output(figure_id, "figure"),
            Output(f"{figure_id}-error", "is_open"),
            Output(f"{figure_id}-error-message", "children"),
        ]
        for figure_id in sorted(PLOTS.keys())
    ]

    @app.callback(
        [
            Output("toast-error", "is_open"),
            Output("toast-error-message", "children"),
            *[out for figure in figure_outputs for out in figure],
        ],
        [Input("submit", "n_clicks")],
    )
    def draw_plots(submit):
        """
        Creates the needed plots using functions in dynophores.webapp.plots.

        Parameters
        ----------
        TODO
            TODO

        Returns
        -------
        tuple of length TODO
            bool: True if error message should be displayed, False otherwise
            str or None: Error message, if any
            For each figure, append:
                dict: figure data for plotly Graph
                bool: True if error message should be displayed for this figure, False otherwise
                str or None: Error message for this figure, if any
        """

        if submit is None:
            # Return is open/closed, error_msg
            return [False, None] + [{}, False, None] * len(PLOTS)

        # First, get all data
        try:
            dynophore = load_dynophore()
        except Exception:  # if error occured, print in general error box
            logger.error(traceback.format_exc())
            error_msg = "\n".join(traceback.format_exc(limit=0, chain=False).splitlines()[1:])
            return [True, error_msg] + [{}, False, None] * len(PLOTS)

        # Now, try to plot figures
        result = [False, None]
        for figure_id, meta in sorted(PLOTS.items()):
            try:
                figure = meta["plotter"](
                    dynophore, title=meta["title"], yaxis=meta["yaxis"],
                )
                error = False
                error_msg = None
            except Exception:
                logger.error(traceback.format_exc())
                figure = {}
                error = True
                error_msg = "\n".join(traceback.format_exc(limit=0, chain=False).splitlines()[1:])
            finally:
                result.extend([figure, error, error_msg])

        return result

    @app.callback(
        Output('output-data-upload', 'children'),
        [
            Input('upload-data', 'contents')
        ],
        [
            State('upload-data', 'filename')
        ]
    )
    def update_output(list_of_contents, list_of_filenames):
        """
        Update output from data upload.

        Parameters
        ----------
        list_of_contents : list of str
            List of file contents.
        list_of_filenames : list of str
            List of filenames.

        Returns
        -------
        list of dash_html_components.Div.Div
            List of parsed file contents.
        """

        if list_of_contents is not None:

            children = [
                parse_contents(c, n) for c, n in
                zip(list_of_contents, list_of_filenames)]

            return children

    @app.callback(
        [
            Output("dynophore-3d-collapse", "is_open"),
            Output("dynophore-2d-collapse", "is_open"),
            Output("interaction-heatmap-collapse", "is_open"),
            Output("superfeatures-collapse", "is_open"),
            Output("interaction-partners-collapse", "is_open"),
        ],
        [
            Input("dynophore-3d-toggle", "n_clicks"),
            Input("dynophore-2d-toggle", "n_clicks"),
            Input("interaction-heatmap-toggle", "n_clicks"),
            Input("superfeatures-toggle", "n_clicks"),
            Input("interaction-partners-toggle", "n_clicks"),
        ],
        [
            State("dynophore-3d-collapse", "is_open"),
            State("dynophore-2d-collapse", "is_open"),
            State("interaction-heatmap-collapse", "is_open"),
            State("superfeatures-collapse", "is_open"),
            State("interaction-partners-collapse", "is_open"),
        ],
    )
    def toggle_collapse(n1, n2, n3, n4, n5, is_open1, is_open2, is_open3, is_open4, is_open5):
        """
        Open/close cards upon toggle.

        Parameters
        ----------
        n1 : int
            Number of clicks.
        is_open1 : bool
            Card is open? (input)

        Returns
        -------
        list of bool
            Card is open? (output)
        """

        ctx = dash.callback_context
        if not ctx.triggered:
            return False, False, False, False, False
        else:
            button_id = ctx.triggered[0]["prop_id"].split(".")[0]

        # Return state (open=True, closed=False)
        if button_id == "dynophore-3d-toggle" and n1:
            return not is_open1, is_open2, is_open3, is_open4, is_open5
        elif button_id == "dynophore-2d-toggle" and n2:
            return is_open1, not is_open2, is_open3, is_open4, is_open5
        elif button_id == "interaction-heatmap-toggle" and n3:
            return is_open1, is_open2, not is_open3, is_open4, is_open5
        elif button_id == "superfeatures-toggle" and n4:
            return is_open1, is_open2, is_open3, not is_open4, is_open5
        elif button_id == "interaction-partners-toggle" and n5:
            return is_open1, is_open2, is_open3, is_open4, not is_open5
        return True, True, True, True, True

    return app
