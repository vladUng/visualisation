import dash_core_components as dcc
import dash_html_components as html
from dash.dependencies import Input, Output, State
import dash_bio as dashbio


from plotlyflask.plotlydash.main import menu
from plotlyflask.plotlydash.main import app as dash_app

import numpy as np
import pandas as pd

from os import path, walk
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

def create_dropdown(files, plot_type):

    return html.Div(children=[
        html.H6('Select a file'),
        dcc.Dropdown(
                id="select-file-{}".format(plot_type), value=files[0],
                options=[{"label": i, "value": i} for i in files]
        ),
    ])

def create_config_vulcano(files):

    return html.Div(id="config-vulcano", children=[
        html.H4('Fold Changes sizes'),
        create_dropdown(files, "volcano"),
        html.Br(),
        dcc.RangeSlider(
            id='default-volcanoplot-input',
            min=-3,
            max=3,
            step=0.05,
            marks={i: {'label': str(i)} for i in range(-3, 3)},
            value=[-1.0, 1]
        ),
        html.Br(),
        html.Button(id='plot-volcano', n_clicks=0, children='Plot Volcano'),
        html.Br(),
    ])

def create_config_piplot(files):

    return html.Div(id="config-pi", style={"height": "150pt"}, children=[
        html.H4('Settings for π plot'),
        create_dropdown(files, "pi"),
        html.Br(),
        html.Button(id='plot-pi', n_clicks=0, children='Plot π plot'),
    ])

def init_callbacks(dash_app):
    @dash_app.callback(
        [Output("volcano-text-output", "children"),
         Output('figure-volcano', "figure")],
        [Input('plot-volcano', 'n_clicks'),
         Input('default-volcanoplot-input', 'value')],
        [State("select-file-volcano", "value")]
    )
    def plotVolcano(btn, fold_changes, filename):
        ret_string = ""
        data_dict = import_data("data/VolcanoPlots/" + filename)

        figure = draw_volcano(data_dict["data"], fold_changes)
        return ret_string, figure

    @dash_app.callback(
    [Output("pi-plot-text-output", "children"),
        Output('figure-pi-plot', "figure")],
    [Input('plot-pi', 'n_clicks')],
    [State("select-file-pi", "value")]
    )
    def plotPi(btn, filename):
        ret_string = ""
        data_dict = import_data("data/VolcanoPlots/" + filename)

        figure = draw_volcano(data_dict["data"], [-1, 1])
        return ret_string, figure



filename = "Basal_split.tsv"
#filename = "BasalSplit_vsOthers_v2_vulcano.tsv"
filename = "BasalSplit_v2_vulcano.tsv"

files = next(walk("data/VolcanoPlots/"), (None, None, []))[2]

init_callbacks(dash_app)

# Create Layout
layout = html.Div(children=[
    menu,
    html.H1(children='Volcano plot'),
    html.Div(id="parrent-all-volcano", children=[
        html.Hr(),
        html.Div(id="figure-contols-volcano",  style={"column-count": "2"}, children=[
            create_config_vulcano(files),
            create_config_piplot(files)
        ]),
        html.Hr(),
        html.Br(),
        html.Div(id="figures",  style={"column-count": "2"}, children=[
            html.Div(id='volcano-text-output'),
            dcc.Graph(
                id='figure-volcano',
                figure={"data": {},
                        "layout": {"title": "No Volcano displayed", "height": 800}
                        }),
            html.Div(id='pi-plot-text-output'),
            dcc.Graph(
                id='figure-pi-plot',
                figure={"data": {},
                        "layout": {"title": "No π plot displayed", "height": 800}
                        }),            
        ])
    ]),
])
