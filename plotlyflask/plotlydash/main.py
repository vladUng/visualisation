import dash
from flask import Flask
from dash import dcc
from dash import html

"""Construct core Flask application."""
server = Flask(__name__, instance_relative_config=True)

external_stylesheets = ['https://codepen.io/chriddyp/pen/bWLwgP.css']
app = dash.Dash(
    server=server,
    external_stylesheets=external_stylesheets,
    suppress_callback_exceptions=True
)

menu = html.Div([
    html.Br(),
    dcc.Link('Go to Gene Visualisation', href='/gene-vis'), html.Br(),
    dcc.Link('Go to Gene Differentiation', href='/gene-diff'), html.Br(),
    dcc.Link('Go to Umap', href='/manyfold'), html.Br(),
    dcc.Loading(children=[html.Div(id="ls-loading-output-1")], type="default"),
    html.Br()
])