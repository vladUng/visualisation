import time
from os import path

import numpy as np
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
import umap
from dash import dcc, html
from dash.dependencies import Input, Output, State

from plotlyflask.plotlydash.main import app as dash_app
from plotlyflask.plotlydash.main import menu

pd.options.plotting.backend = "plotly"

def import_data(base_path, filename="baker_BK.tsv", meta_filename="meta_baker_BK.tsv"):
    """Creates the dataframe used for the tool and the metadata

    Args:
        base_path ([String]): Path to the data

    Returns: 
        DataDict: A dictionary which contains the the all_tsv and metadata 
    """
    start_time = time.time()
    if path.exists(base_path):

        ret_dict = {}
        ret_dict["data"] = pd.read_csv(base_path + filename, delimiter="\t")
        ret_dict["data"] = ret_dict["data"].set_index("genes")
        ret_dict["metadata"] = pd.read_csv(base_path + meta_filename, delimiter="\t")
        print("Finished loading the data in {}".format(time.time()-start_time))
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

    for cluster_size in range(4, 6):
        df["RawKMeans_CS_{}".format(cluster_size)] = df_meta["RawKMeans_CS_{}".format(
            cluster_size)].astype(str)
        df["Ward_CS_{}".format(cluster_size)] = df_meta["Ward_CS_{}".format(
            cluster_size)].astype(str)
        df["GaussianMixture_CS_{}".format(
            cluster_size)] = df_meta["GaussianMixture_CS_{}".format(cluster_size)].astype(str)
        df["FuzzyCMeans_CS_{}".format(
            cluster_size)] = df_meta["GaussianMixture_CS_{}".format(cluster_size)].astype(str)


def draw_umap(data, meta_df, n_neighbors=15, min_dist=0.1, n_components=2,metric='cosine', title='', colour="TCGA408_classifier"):
    umap_model = umap.UMAP(n_neighbors=n_neighbors, min_dist=min_dist, n_components=n_components, metric=metric, )

    u_fit = umap_model.fit_transform(data.values)
    columns = []
    if n_components == 3:
        columns = ["UMAP 1", "UMAP 2", "UMAP 3"]
    else:
        columns = ["UMAP 1", "UMAP 2"]

    umap_df = pd.DataFrame(u_fit, columns=columns)
    dummy_meta = meta_df[meta_df["Samples"].isin(data.index.values)].rename(columns={"Samples": "sample"}).sort_values(by="sample").reset_index(drop=True)

    tcga_add_metadata(umap_df, dummy_meta)

    layout = go.Layout(title="UMAP plot TPM values", height=700)

    fig = {}
    if n_components == 3:
        fig = px.scatter_3d(umap_df, x="UMAP 1", y="UMAP 2", z="UMAP 3", color=colour, title=title)
    else:

        fig = px.scatter(umap_df, x="UMAP 1", y="UMAP 2", color=colour, title=title)

    fig.update_layout(layout)

    return fig


def init_callbacks(dash_app, data_dict):

    # For Figure 1
    @dash_app.callback(
        [Output('umap-config-1', "children"), Output("umap-plot-1", "figure")],
        [Input("compute-viz-1", 'n_clicks')],
        [State("neighbours-slider-1", 'value'),
         State("distance-slider-1", 'value'), State("select-metrics-1","value"), State("select-colouring-1", "value"),
         State("select-components-1", "value")]
    )
    def plotUMAP(btn, neighbours, distance, metric, colouring, n_comp):

        ret_string = "Umap with the following config: Neighbours: {}, Distance: {}, Metric: {}, Colour: {}, Components: {}".format(
            neighbours, distance, metric, colouring, n_comp)

        figure = {}
        if btn:
            figure = draw_umap(data_dict["data"], data_dict["metadata"], n_neighbors=neighbours, min_dist=distance, n_components=n_comp, metric=metric, colour=str(colouring))

        return ret_string, figure

    @dash_app.callback(
        Output('neighbours-output-container-1', "children"),
        [Input("neighbours-slider-1", 'value')]
    )
    def updateNeighboursSlider(value):
        return "Selected neighbours value {}".format(value)

    @dash_app.callback(
        Output('distance-output-container-1', "children"),
        [Input("distance-slider-1", 'value')]
    )
    def updateDistanceSlider(value):
        return "Selected distance value {}".format(value)

    # For Figure 2
    @dash_app.callback(
        [Output('umap-config-2', "children"), Output("umap-plot-2", "figure")],
        [Input("compute-viz-2", 'n_clicks')],
        [State("neighbours-slider-2", 'value'),
         State("distance-slider-2", 'value'), State("select-metrics-2", "value"), State("select-colouring-2", "value"),
         State("select-components-2", "value")]
    )
    def plotUMAP_2(btn, neighbours, distance, metric, colouring, n_comp):

        ret_string = "Umap with the following config: Neighbours: {}, Distance: {}, Metric: {}, Colour: {}, Components: {}".format(
            neighbours, distance, metric, colouring, n_comp)

        figure = {}
        if btn:
            figure = draw_umap(data_dict["data"], data_dict["metadata"], n_neighbors=neighbours,
                            min_dist=distance, n_components=n_comp, metric=metric, colour=str(colouring))
        return ret_string, figure

    @dash_app.callback(
        Output('neighbours-output-container-2', "children"),
        [Input("neighbours-slider-2", 'value')]
    )
    def updateNeighboursSlider(value):
        return "Selected neighbours value {}".format(value)

    @dash_app.callback(
        Output('distance-output-container-2', "children"),
        [Input("distance-slider-2", 'value')]
    )
    def updateDistanceSlider(value):
        return "Selected distance value {}".format(value)


def create_config_menu(no_figure, df_meta):
    default_comp = 2
    default_metric = "cosine"

    if no_figure == 2:
        default_comp = 3
        default_metric = "euclidean"

    colouring_options = ["TCGA408_classifier", "2019_consensus_classifier", "TCGA_2017_AM_remap"]
    colouring_options = colouring_options + list(df_meta.columns[-8:])

    return html.Div(id="config-{}".format(no_figure),  # style={"width": "50%", "column-count": "2"},
                    children=[
        html.H5("Choose the parameters for Figure {}".format(no_figure)),
        html.H6("Neighbours Slider"),
        dcc.Slider(
            id="neighbours-slider-{}".format(no_figure),
            min=0, max=100, step=5,
            marks={0: '0', 25: '25', 50: '50', 75: '75', 100: '100'},
            value=10
        ),
        html.Div(id='neighbours-output-container-{}'.format(no_figure)),
        html.H6("Distance Slider"),
        dcc.Slider(
            id="distance-slider-{}".format(no_figure),
            min=0, max=1, step=0.1,
            marks={0: '0', 0.25: '0.25', 0.5: '0.5',
                   0.75: '0.75', 1.0: '1.0'},
            value=0.2
        ),
        html.Div(id="distance-output-container-{}".format(no_figure)),
        html.Br(),
        html.H6('Select number of components'),
        dcc.RadioItems(
            id="select-components-{}".format(no_figure), value=default_comp,
            options=[
                {'label': 'Two', 'value': 2, },
                {'label': 'Three', 'value': 3},
            ]
        ),
        html.H6('Select colour'),
        dcc.Dropdown(
            id="select-colouring-{}".format(no_figure), value='TCGA_2017_AM_remap',
            options=[{"label": i, "value": i} for i in colouring_options]
        ),
        html.H6('Select distance metric'),
        dcc.RadioItems(
            id="select-metrics-{}".format(no_figure), value=default_metric,
            options=[
                {'label': 'Cosine', 'value': 'cosine', },
                {'label': 'Euclidean', 'value': 'euclidean'},
                {'label': 'Manhattan', 'value': 'manhattan'},
            ]
        ),
        html.Br(),
        html.Button(id="compute-viz-{}".format(no_figure),
                    n_clicks=0, children='Figure {}'.format(no_figure)),
        html.Br(),
    ])

# data_dict = import_data("JBU_data/UMAP/")
data_dict = import_data("data/gene_viz/")
init_callbacks(dash_app, data_dict)


# Create Layout
layout = html.Div(children=[
    menu,
    html.H1(children='UMAP visualisation tool'),
    html.Div(id="parrent-all-umap", children=[
        html.Hr(),
        html.Div(id="figure-contols-umap", style={"column-count": "2"}, children=[
            create_config_menu(1, data_dict["metadata"]),
            create_config_menu(2, data_dict["metadata"]),
        ]),
        html.Hr(),
        html.Div(id="figures",  style={"column-count": "2"}, children=[
            html.Div(id='umap-config-{}'.format(1)),
            dcc.Graph(
                id='umap-plot-{}'.format(1),
                figure={"data": {},
                        "layout": {"title": "Select the configuration first",
                                   "height": 300}
                        }),
            html.Div(id='umap-config-{}'.format(2)),
            dcc.Graph(
                id='umap-plot-{}'.format(2),
                figure={"data": {},
                        "layout": {"title": "Select the configuration first",
                                   "height": 300}
                        })
        ])
    ]),
])
