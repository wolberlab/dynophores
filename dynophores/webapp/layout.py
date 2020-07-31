"""
This module defines the layout and style of the webapp.
"""

from functools import partial

import dash_html_components as html
import dash_core_components as dcc
import dash_bootstrap_components as dbc

from .plots import boxplot

THEME = dbc.themes.JOURNAL


def layout():
    """
    Concatenates all component groups as defined in the functions above. This is
    the final result that will be rendered in the webapp.
    """

    return dbc.Container(
        [
            *_header(),
            *_upload(),
            *_submit_button(),
            *_general_errors(),
            *_interaction_heatmap(),
            *_superfeatures(),
            *_interaction_partners(),
        ]
    )


def _header():
    """
    Creates the top part of the webapp.
    """
    return (
        dbc.Alert(
            "The app is WIP.", color="warning"
        ),
        html.H1("Dynophore Viewer"),
        html.Div("Analyze dynamic protein-ligand interactions"),
        html.Hr(),
    )


def _upload():
    """
    Create an upload button.

    Components
    ----------
    * `upload-data`: Dataset upload.
    * `output-data-upload`: Output from dataset upload.
    """

    return (
        dbc.Row(
            [
                dbc.Col([dbc.Label("Choose a dataset:")], width=2),
                dbc.Col(
                    [
                        dcc.Upload(
                            id='upload-data',
                            children=html.Div([
                                'Drag and Drop or ',
                                html.A('Select Files')
                            ]),
                            style={
                                'width': '100%',
                                'height': '60px',
                                'lineHeight': '60px',
                                'borderWidth': '1px',
                                'borderStyle': 'dashed',
                                'borderRadius': '5px',
                                'textAlign': 'center',
                                'margin': '10px'
                            },
                            # Allow multiple files to be uploaded
                            multiple=True
                        ),
                        html.Div(id='output-data-upload'),
                    ]
                ),
            ]
        ),
    )


def _submit_button():
    """
    Create submit button.

    Components
    ----------
    * `submit`: Submit button.
    """

    return (
        [
            dbc.Button(
                id='submit',
                children="Submit",
                color="primary",
                style={'margin': '10px'}
            )
        ]
    )


def _general_errors():
    """
    Toast with an error message.

    Components
    ----------
    * `toast-error`
    * `toast-error-message`
    """
    return (
        dbc.Toast(
            [html.P(id="toast-error-message")],
            id="toast-error",
            header="An error occured!",
            icon="danger",
            dismissable=True,
            is_open=False,
            style={"max-width": "100%"},
        ),
    )


def _interaction_heatmap():
    """
    Creates the interaction heatmap area for the dynophore.

    Components
    ----------
    * `interaction-heatmap`: plotly graph for interaction heatmap (+ their errors; check `_plot_area`)
    * `interaction-heatmap-toggle`: button controlling open/closed state of plot area
    * `interaction-heatmap-collapse`: collapsible area that contains plots
    """

    return (
        dbc.Card(
            [
                dbc.CardHeader(
                    html.H2(
                        dbc.Button(
                            "Interaction heatmap", color="link", id="interaction-heatmap-toggle",
                        )
                    )
                ),
                dbc.Collapse(
                    [*_plot_area("interaction-heatmap")],
                    id="interaction-heatmap-collapse",
                ),
            ]
        ),
    )


def _superfeatures():
    """
    Creates the superfeatures area for the dynophore.

    Components
    ----------
    * `superfeatures`: plotly graph for superfeatures occurrences (+ their errors; check `_plot_area`)
    * `superfeatures-toggle`: button controlling open/closed state of plot area
    * `superfeatures-collapse`: collapsible area that contains plots
    """

    return (
        dbc.Card(
            [
                dbc.CardHeader(
                    html.H2(
                        dbc.Button(
                            "Superfeatures", color="link", id="superfeatures-toggle",
                        )
                    )
                ),
                dbc.Collapse(
                    [*_plot_area("superfeatures")],
                    id="superfeatures-collapse",
                ),
            ]
        ),
    )


def _interaction_partners():
    """
    Creates the superfeatures area for the dynophore.

    Components
    ----------
    * `interaction-partners`: plotly graph for interaction partner analysis (+ their errors; check `_plot_area`)
    * `interaction-partners-toggle`: button controlling open/closed state of plot area
    * `interaction-partners-collapse`: collapsible area that contains plots
    """

    return (
        dbc.Card(
            [
                dbc.CardHeader(
                    html.H2(
                        dbc.Button(
                            "Interaction partners", color="link", id="interaction-partners-toggle",
                        )
                    )
                ),
                dbc.Collapse(
                    [*_plot_area("interaction-partners")],
                    id="interaction-partners-collapse",
                ),
            ]
        ),
    )


def _plot_area(identifier):
    """
    Loading area for plot and some errors if needed.

    Components
    ----------
    * `graph-{identifier}-error`: Box containing error info, if needed
    * `graph-{identifier}-error-message`: Error message
    * `graph-{identifier}-loading`: Loading layer
    * `{identifier}`: The target dcc.Graph component
    """
    return (
        dbc.Toast(
            [html.P(id=f"{identifier}-error-message")],
            id=f"{identifier}-error",
            header="An error occured!",
            icon="danger",
            dismissable=True,
            is_open=False,
            style={"max-width": "100%"},
        ),
        dcc.Loading(
            id=f"{identifier}-loading", children=[dcc.Graph(id=identifier)], type="default",
        ),
    )


PLOTS = {
    "interaction-heatmap": {
        "title": "Interaction occurrences: Superfeature vs. interaction partners",
        "yaxis": "TBA",
        "plotter": partial(boxplot, selector=lambda dynophore: dynophore.count),
    },
    "superfeatures": {
        "title": "Superfeature occurrences",
        "yaxis": "TBA",
        "plotter": partial(boxplot, selector=lambda dynophore: dynophore.superfeatures_occurrences),
    },
    "interaction-partners": {
        "title": "Interaction partner occurrences and distances",
        "yaxis": "TBA",
        "plotter": partial(boxplot, selector=lambda dynophore: dynophore.envpartner_occurrences),
    },
}
