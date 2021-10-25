import dash
from flask import Flask
import dash_html_components as html
import dash_core_components as dcc

"""Construct core Flask application."""
server = Flask(__name__, instance_relative_config=True)

menu = html.Div([
    html.Br(),
    dcc.Link('Go to Gene Visualisation', href='/gene-vis'), html.Br(),
    dcc.Link('Go to Gene Differentiation', href='/gene-diff'), html.Br(),
    dcc.Link('Go to Umap', href='/manyfold'), html.Br(),
    dcc.Loading(children=[html.Div(id="ls-loading-output-1")], type="default"),
    html.Br()
])

external_stylesheets = ['https://codepen.io/chriddyp/pen/bWLwgP.css']
app = dash.Dash(
    server=server,
    # routes_pathname_prefix='/',
    external_stylesheets=external_stylesheets,
    suppress_callback_exceptions=True
)
