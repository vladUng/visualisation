import dash_core_components as dcc
import dash_html_components as html
from dash.dependencies import Input, Output, State
import dash_bio as dashbio
import plotly.express as px


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

def draw_pi_plot(df_1, df_2, file_1, file_2):

    title = ""
    # Decide which is bigger, the number of genes may differ
    if df_1.shape[0] < df_2.shape[0]:
        first_df = df_2.copy(deep=True).set_index("genes")
        second_df = df_1.copy(deep=True).set_index("genes")
        title = "X: {} on  vs Y: {}".format(file_2, file_1)
    else:
        first_df = df_1.copy(deep=True).set_index("genes")
        second_df = df_2.copy(deep=True).set_index("genes")
        title = "X: {} on  vs Y: {}".format(file_1, file_2)
        
    # compute the values
    first_df["x"] = -np.log10(first_df["q"]) * first_df["fold_change"]
    second_df["y"] = -np.log10(second_df["q"]) * second_df["fold_change"]

    # Note the genes number may differ and we're setting the index of the DataFrame which has the most genes (i.e. rows)
    #  However, there might be some genes in the second df which are not in the first one. Regardless, we set the nan values to 0 (points will be added to the center)
    dummy_df = pd.concat([first_df["x"], second_df["y"]], axis=1).fillna(0).reset_index().rename(columns={"index":"genes"})

    fig = px.scatter(dummy_df, x="x", y="y", hover_data=["genes"], title=title)


    fig.add_shape(type='line',
                    x0=dummy_df["x"].min()*1.5, y0=0,
                    x1=dummy_df["x"].max()*1.5, y1=0,
                    line=dict(color='Black', width=1),
                    xref='x', yref='y')

    fig.add_shape(type='line',
                    x0=0, y0=dummy_df["y"].min()*1.5,
                    x1=0, y1=dummy_df["y"].max()*1.5,
                    line=dict(color='Black', width=1),
                    xref='x',
                    yref='y')

    # fig.update_traces(marker=dict(size=12, line=dict(width=2, color='DarkSlateGrey')), selector=dict(mode='markers'))

    return fig

def create_dropdown(files, plot_type, text="Select a file"):

    return html.Div(children=[
        html.H6(text),
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
        create_dropdown(files, "pi-1", text="π DEA 1. Select a file."),
        create_dropdown(files, "pi-2", text="π DEA 2. Select a file."),
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
    [State("select-file-pi-1", "value"), State("select-file-pi-2", "value")]
    )
    def plotPi(btn, file_1, file_2):
        ret_string = ""
        data_dict_1 = import_data("data/VolcanoPlots/" + file_1)
        data_dict_2 = import_data("data/VolcanoPlots/" + file_2)

        figure = draw_pi_plot(data_dict_1["data"], data_dict_2["data"], file_1, file_2)
        return ret_string, figure



files = next(walk("data/VolcanoPlots/"), (None, None, []))[2]

init_callbacks(dash_app)

# Create Layout
layout = html.Div(children=[
    menu,
    html.H1(children='Volcano plot'),
    html.Div(id="parrent-all-volcano", children=[
        html.Hr(),
        html.Div(id="figure-contols-volcano",  style={"column-count": "1"}, children=[
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
