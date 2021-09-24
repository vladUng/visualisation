import pandas as pd
import dash
import dash_html_components as html
import dash_core_components as dcc
from dash.dependencies import Input, Output, State
import plotly.express as px
import plotly.graph_objects as go

import numpy as np

from os import path, name
import time
from datetime import datetime

import umap


def import_data(base_path):
    """Creates the dataframe used for the tool and the metadata

    Args:
        base_path ([String]): Path to the data

    Returns: 
        DataDict: A dictionary which contains the the all_tsv and metadata 
    """
    start_time = time.time()
    if path.exists(base_path):

        ret_dict = {}
        ret_dict["data"] = pd.read_csv(base_path + "TPM.tsv", delimiter="\t")
        ret_dict["data"] = ret_dict["data"].set_index("Unnamed: 0")
        ret_dict["metadata"] = pd.read_csv( base_path + "metadata.tsv", delimiter="\t")

        # TPM and metadata
        return ret_dict
    else:
        return None


def tcga_add_metadata(df, df_meta):

    # first numeric data
    for col in ["cigarettes_per_day", "bmi", "weight", "height"]:
        df_meta[col] = df_meta[col].replace("--", np.nan).astype(np.float)
        df_meta[col].fillna(df_meta[col].median(), inplace=True)
        df[col] = df_meta[col].values

    df["2019_consensus_classifier"] = df_meta["2019_consensus_classifier"].values
    df["gender"] = df_meta["gender"].values
    df["TCGA408_classifier"] = df_meta["TCGA408_classifier"].values
    df["samples"] = df_meta["sample"].values

    df["2019_consensus_classifier"] = df_meta["2019_consensus_classifier"]
    df["TCGA_2017_AM_remap"] = df_meta["TCGA_2017_AM_remap"]


def draw_umap(data, meta_df, n_neighbors=15, min_dist=0.1, n_components=2,
              metric='cosine', title='', colour="TCGA408_classifier"):
    umap_model = umap.UMAP( n_neighbors=n_neighbors, min_dist=min_dist, n_components=n_components, metric=metric)

    u_fit = umap_model.fit_transform(data.values)
    umap_df = pd.DataFrame(u_fit, columns=["UMAP 1", "UMAP 2"])

    dummy_meta = meta_df[meta_df["Samples"].isin(data.index.values)].rename(columns={"Samples": "sample"}).sort_values(by="sample").reset_index(drop=True)

    tcga_add_metadata(umap_df, dummy_meta)

    layout = go.Layout(title="UMAP plot TPM values", height=700)

    fig = px.scatter(umap_df, x="UMAP 1", y="UMAP 2", color=colour, title=title)
    fig.update_layout(layout)

    return fig


def init_callbacks(dash_app, data_dict):

    @dash_app.callback(
        [Output('umap-config', "children"), Output("umap-plot", "figure")],
        [Input("compute-viz", 'n_clicks')],
        [State("neighbours-slider", 'value'),
         State("distance-slider", 'value'), State("select-metrics", "value"), State("select-colouring", "value")]
    )
    def plotUMAP(btn, neighbours, distance, metric, colouring):

        n_comp = 2
        # metric = 'cosine'
        colouring = "TCGA408_classifier"
        ret_string = "Umap with the following config: neighbours: {}, Distance: {}, Metric: {}, Colour: {}, Components: {}".format(
            neighbours, distance, metric, colouring, n_comp)

        figure = draw_umap(data_dict["data"], data_dict["metadata"], n_neighbors=neighbours,
                           min_dist=distance, n_components=n_comp, metric=metric, colour=str(colouring))

        return ret_string, figure

    @dash_app.callback(
        Output('neighbours-output-container', "children"),
        [Input("neighbours-slider", 'value')]
    )
    def updateNeighboursSlider(value):
        return "Selected neighbours value {}".format(value)

    @dash_app.callback(
        Output('distance-output-container', "children"),
        [Input("distance-slider", 'value')]
    )
    def updateDistanceSlider(value):
        return "Selected distance value {}".format(value)

def create_config_menu():
    return html.Div(id="config", style={"width": "50%", "column-count": "2"},
              children= [
                html.H6("Neighbours Slider"),
                dcc.Slider(
                    id='neighbours-slider',
                    min=0, max=100, step=10,
                    marks={0: '0', 25: '25', 50: '50', 75: '75', 100: '100'},
                    value=20
                ),
                html.Div(id='neighbours-output-container'),
                html.H6("Distance Slider"),
                dcc.Slider(
                    id='distance-slider',
                    min=0, max=1, step=0.1,
                    marks={0: '0', 0.25: '0.25', 0.5: '0.5', 0.75: '0.75', 1.0: '1.0'},
                    value=0.5
                ),
                html.Div(id='distance-output-container'),
                html.H6('Select colour'),
                dcc.RadioItems(
                    id="select-colouring", value='TCGA_2017_AM_remap',
                    options = [
                        {'label':'TCGA 2017', 'value':'TCGA408_classifier',},
                        # {'label':'Remap TCGA 2017', 'value':'TCGA_2017_AM_remap'},
                        # {'label': 'Remap Consensus', 'value': '2019_consensus_classifier'},
                    ]
                ),
                html.H6('Select distance metric'),
                dcc.RadioItems(
                    id="select-metrics", value='cosine',
                    options = [
                        {'label':'Cosine', 'value':'cosine',},
                        {'label':'Euclidean', 'value':'euclidean'},
                        {'label': 'Manhattan', 'value': 'manhattan'},
                    ]
                ),
            ])

def init_dashboard(server):
    """Create a Plotly Dash dashboard."""
    external_stylesheets = ['https://codepen.io/chriddyp/pen/bWLwgP.css']
    dash_app = dash.Dash(
        server=server,
        routes_pathname_prefix='/manyfold/',
        external_stylesheets=external_stylesheets
    )

    data_dict = import_data("data/UMAP/")
    init_callbacks(dash_app, data_dict)

    # Create Layout
    dash_app.layout = html.Div(children=[
        html.H1(children='UMAP visualisation tool'),
        html.Div(id="parrent-all", children=[
            html.Hr(),
            html.H5("Choose the parameters"),
            create_config_menu(),
            html.Br(),
            html.Button(id='compute-viz', n_clicks=0, children='Compute'),
            html.Br(),
            html.Hr(),
            html.Div(id='umap-config'),
            dcc.Graph(
                id='umap-plot',
                figure={ "data": {}, 
                "layout": {"title": "Select the configuration first", 
                "height": 300} })
        ]),
    ])
    return dash_app.server
