import dash_core_components as dcc
import dash_html_components as html
from dash.dependencies import Input, Output

from plotlyflask.plotlydash.main import app

menu = html.Div([
    html.Div(id='gene-viz-dislay'),
    dcc.Link('Go to Gene Visualisation', href='/gene-vis'),
    html.Div(id='gene-diff-dislay'),
    dcc.Link('Go to Gene Differentiation', href='/gene-diff'),
    html.Div(id='manyfold'),
    dcc.Link('Go to Umap', href='/manyfold'),
])


layout = html.Div([
    menu,
    "TBD"
])
