# from .features import gene_diff as gf
# from .features import manyfold as mf
# from .features import gene_viz as gv

import dash
from flask import Flask

# from plotlyflask.plotlydash.features import gene_diff as gf
# from plotlyflask.plotlydash.features import gene_viz as gv
# from plotlyflask.plotlydash.features import manyfold as mf


"""Construct core Flask application."""
server = Flask(__name__, instance_relative_config=True)

external_stylesheets = ['https://codepen.io/chriddyp/pen/bWLwgP.css']
app = dash.Dash(
    server=server,
    # routes_pathname_prefix='/',
    external_stylesheets=external_stylesheets
)
