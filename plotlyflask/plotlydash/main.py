import dash
from flask import Flask

"""Construct core Flask application."""
server = Flask(__name__, instance_relative_config=True)

external_stylesheets = ['https://codepen.io/chriddyp/pen/bWLwgP.css']
app = dash.Dash(
    server=server,
    # routes_pathname_prefix='/',
    external_stylesheets=external_stylesheets
)
