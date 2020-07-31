"""
This module defines the webapp, the Dynophore Viewer.
"""

from dash import Dash

from .layout import layout, THEME
from .callbacks import register_callbacks


def create_app():
    """
    Bootstrap Dash app with layout, theme and callbacks.
    """

    app = Dash(__name__, external_stylesheets=[THEME])
    app.title = "Dynophores"
    app.layout = layout()
    register_callbacks(app)
    return app


def main():
    app = create_app()
    app.run_server(debug=True)


if __name__ == "__main__":
    main()
