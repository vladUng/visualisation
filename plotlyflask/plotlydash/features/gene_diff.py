import dash_core_components as dcc
import dash_html_components as html
from dash.dependencies import Input, Output
import dash_bio as dashbio


from plotlyflask.plotlydash.main import menu
from plotlyflask.plotlydash.main import app as dash_app

import numpy as np
import pandas as pd

from os import path
import time

##### Functions for .tsv handling #####


def import_data(fullPath):
    start_time = time.time()
    if path.exists(fullPath):
        ret_dict = {}
        ret_dict["data"] = pd.read_csv(fullPath, delimiter="\t")
        print("Finished loading the data in {}".format(time.time()-start_time))
        return ret_dict
    else:
        return None


def draw_volcano(df, fold_changes):
    return dashbio.VolcanoPlot(
        dataframe=df,
        effect_size="fold_change",
        gene="genes",
        snp=None,
        p="p",
        genomewideline_value=2.5,
        logp=True,
        effect_size_line=fold_changes,
        ylabel="-log10(q)",
        xlabel="log2(FC)",
        col='#2A3F5F',
        point_size=8,
        effect_size_line_width=4,
        genomewideline_width=2
    )


def init_callbacks(dash_app, data_dict):
    @dash_app.callback(
        [Output("volcano-text-output", "children"),
         Output('figure-volcano', "figure")],
        [Input('plot-volcano', 'n_clicks'),
         Input('default-volcanoplot-input', 'value')]
    )
    def plotVolcano(btn, fold_changes):
        ret_string = ""
        figure = draw_volcano(data_dict["data"], fold_changes)
        return ret_string, figure


filename = "Basal_split.tsv"
data_dict = import_data("data/VolcanoPlots/" + filename)
init_callbacks(dash_app, data_dict)

# Create Layout
layout = html.Div(children=[
    menu,
    html.H1(children='Volcano plot'),
    html.Div(id="parrent-all-volcano", children=[
        html.Hr(),
        html.Div(id="figure-contols-volcano", children=[
            'Fold Changes sizes',
            dcc.RangeSlider(
                id='default-volcanoplot-input',
                min=-3,
                max=3,
                step=0.05,
                marks={i: {'label': str(i)} for i in range(-3, 3)},
                value=[-1.0, 1]
            ),
        ]),
        html.Hr(),
        html.Br(),
        html.Div(id='volcano-text-output'),
        dcc.Graph(
            id='figure-volcano',
            figure={"data": {},
                    "layout": {"title": "No Volcano displayed", "height": 800}
                    }),
        html.Button(id='plot-volcano', n_clicks=0, children='Plot'),
    ]),
])
