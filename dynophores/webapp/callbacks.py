"""
This module defines the webapp callbacks.
"""

import base64
import datetime
import io
import logging
import traceback

from dash.dependencies import Input, Output, State
import dash_html_components as html
import dash_table
import pandas as pd

from .layout import PLOTS
from .content import load_dynophore

logger = logging.getLogger(__name__)


def parse_contents(contents, filename, date):
    """
    Parse contents from file upload.

    Parameters
    ----------
    contents : TODO
        TODO
    filename : str
        File name.
    date : str
        Upload date.

    Returns
    -------
    TODO
        TODO
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
    except Exception as e:
        print(e)
        return html.Div([
            'There was an error processing this file.'
        ])

    html_result = html.Div([
        html.H5(filename),
        html.H6(datetime.datetime.fromtimestamp(date)),

        dash_table.DataTable(
            data=df.to_dict('records'),
            columns=[{'name': i, 'id': i} for i in df.columns]
        ),
    ])

    print(type(html_result))

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
            State('upload-data', 'filename'),
            State('upload-data', 'last_modified')
        ]
    )
    def update_output(list_of_contents, list_of_filenames, list_of_dates):
        """
        Update output from data upload.

        Parameters
        ----------
        list_of_contents : TODO
            TODO
        list_of_filenames : list of str
            List of filenames.
        list_of_dates : list of str
            List of dates (last modified).

        Returns
        -------
        TODO
            TODO
        """

        if list_of_contents is not None:

            children = [
                parse_contents(c, n, d) for c, n, d in
                zip(list_of_contents, list_of_filenames, list_of_dates)]

            print(type(children))

            return children
