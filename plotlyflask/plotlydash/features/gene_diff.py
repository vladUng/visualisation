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



styles = {
    'pre': {
        'border': 'thin lightgrey solid',
        'overflowX': 'scroll'
    }
}

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

# Graphs
def draw_volcano(df, fold_changes, selected_data, selected_genes):
    fig = dashbio.VolcanoPlot(
        dataframe=df,
        effect_size="fold_change",
        gene="genes",
        snp=None,
        p="q",
        genomewideline_value=2.5,
        logp=True,
        effect_size_line=fold_changes,
        ylabel="-log10(q)",
        xlabel="log2(FC)",
        col='#2A3F5F',
        point_size=8,
        effect_size_line_width=4,
        genomewideline_width=2,
        annotation="group"
    )
    
    x_col = "fold_change"
    y_col = "-log10(q)"

    # show the selected points
    if selected_data:
        ranges = selected_data['range']
        selection_bounds = {'x0': ranges['x'][0], 'x1': ranges['x'][1],
                            'y0': ranges['y'][0], 'y1': ranges['y'][1]}

        selected_idxs = df[df["genes"].isin(selected_genes)].index

        fig.update_traces(selectedpoints=selected_idxs,
            mode='markers', 
            unselected={'marker': { 'opacity': 0.3 }
            })
    else:
        selection_bounds = {'x0': np.min(df[x_col] -2), 'x1': np.max(df[x_col] + 2),
                            'y0': np.min(df[y_col]-2), 'y1': np.max(df[y_col]) + 2}

    fig = show_selected_genes_vulcano(df, fig)

    fig.update_traces(customdata=df["genes"])
    fig.update_layout(dragmode='select')

    fig.add_shape(dict({'type': 'rect',
                        'line': { 'width': 1, 'dash': 'dot', 'color': 'darkgrey' } },
                       **selection_bounds))
    

    return fig

def show_selected_genes_vulcano(df, fig):
    custom_traces = create_custom_traces()
    colors =  px.colors.qualitative.Vivid
    for idx, trace in enumerate(custom_traces): 
        fig.add_trace(create_gene_trace(df, trace["genes"], name=trace["title"], marker_color=colors[idx]))

    return fig

def show_selected_genes_pi(df_1, df_2, fig):
    custom_traces = create_custom_traces()
    colors =  px.colors.qualitative.Vivid
    colors_2 = px.colors.qualitative.Safe
    for idx, trace in enumerate(custom_traces): 
        fig.add_trace(create_gene_trace(df_1, trace["genes"], name=trace["title"], marker_color=colors[idx], df_2=df_2))

        # fig.add_trace(create_gene_trace(df_2, trace["genes"], name=trace["title"], marker_color=colors[idx], volcano_plot=False, x_axis=False))

    return fig 

def create_custom_traces():

    ifnq_genes = ['B2M', 'BATF2', 'BTN3A3', 'CD74', 'CIITA', 'CXCL10', 'CXCL9',
       'EPSTI1', 'GBP1', 'GBP4', 'GBP5', 'HAPLN3', 'HLA-DPA1', 'IDO1',
       'IFI30', 'IFI35', 'IFIT3', 'IFITM1', 'IL32', 'IRF1', 'NLRC5',
       'PARP9', 'PSMB8-AS1', 'PSMB9', 'SAMHD1', 'SECTM1', 'STAT1', 'TAP1',
       'TNFSF13B', 'TYMP', 'UBE2L6', 'WARS1', 'ZBP1']

    rarg_genes = ['ADGRF1', 'ALDH3B1', 'ANKRD13D', 'ANXA11', 'B3GNT3', 'BHLHE41',
       'CAPN1', 'CDK2AP2', 'CIB1', 'CPTP', 'CTSH', 'GALE', 'HSH2D',
       'MMEL1', 'OAS1', 'PLCD3', 'PNPLA2', 'PRKCD', 'RARG', 'RNPEPL1',
       'TCIRG1', 'THEM6', 'TNFAIP2', 'UBA7', 'UNC93B1', 'VAMP8', 'VSIG2']

    luminal_markers = ["KRT20", "PPARG", "FOXA1", "GATA3", "SNX31", "UPK1A", "UPK2", "FGFR3"]

    basal_markers = ["CD44", "KRT6A", "KRT5", "KRT14", "COL17A1"]

    squamos_markers = ["DSC3", "GSDMC", "TCGM1", "PI3", "TP63"]

    immune_markers = ["CD274", "PDCD1LG2", "IDO1", "CXCL11", "L1CAM", "SAA1"]

    neural_diff = ["MSI1", "PLEKHG4B", "GNG4", "PEG10", "RND2", "APLP1", "SOX2", "TUBB2B"]

    custom_traces = []
    # SB
    custom_traces.append({"genes":ifnq_genes, "title": "SB_IFNQ"})
    custom_traces.append({"genes":rarg_genes, "title": "SB_RARG"})
    # TCGA
    custom_traces.append({"genes":luminal_markers, "title": "TCGA_luminal"})
    custom_traces.append({"genes":basal_markers, "title": "TCGA_basal"})
    custom_traces.append({"genes":squamos_markers, "title": "TCGA_squamos"})
    custom_traces.append({"genes":immune_markers, "title": "TCGA_immune"})
    custom_traces.append({"genes":neural_diff, "title": "TCGA_neuroendocrine"})
    
    return custom_traces

def create_gene_trace(df, genes, name="custom genes", marker_color="yellow", marker_size=8, df_2=None): 

    selected_df = df[df["genes"].isin(genes)]

    x, y = [], []
    if df_2 is None:
        y = -np.log10(selected_df["q"])
        x = selected_df['fold_change']
    else:
        selected_df_2 = df_2[df_2["genes"].isin(genes)]
        
        x = -np.log10(selected_df["q"]) * selected_df["fold_change"]
        y = -np.log10(selected_df_2["q"]) * selected_df_2["fold_change"]

    markers = {"size": marker_size, "color": selected_df.shape[0] * [marker_color] }

    trace = dict(type='scatter', x=x, y= y,  showlegend=True, marker=markers, text=selected_df["genes"], mode="markers", name=name)

    return trace

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
    first_df.rename(columns={"group": "comp_1"}, inplace=True)
    second_df.rename(columns={"group": "comp_2"}, inplace=True)


    dummy_df = pd.concat([first_df[["x", "comp_1"]], second_df[["y", "comp_2"]]], axis=1).fillna(0).reset_index().rename(columns={"index":"genes"})
    
    dummy_df["main"] = "PI_plot"
    fig = px.scatter(dummy_df, x="x", y="y", hover_data=["genes", "comp_1", "comp_2"], color="main", title=title)


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

    fig = show_selected_genes_pi(first_df.reset_index(), second_df.reset_index(), fig)

    return fig

# HTML
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

# Callbacks
def init_callbacks(dash_app):
    @dash_app.callback(
        [Output("volcano-text-output", "children"), Output('click-data', 'children'),
         Output('figure-volcano', "figure")],
        [Input('plot-volcano', 'n_clicks'),
         Input('default-volcanoplot-input', 'value')],  Input('figure-volcano', 'selectedData'),
        [State("select-file-volcano", "value")]
    )
    def plotVolcano(btn, fold_changes, selected_data, filename):
        ret_string = ""
        ret_genes = []
        selected_genes = []

        data_dict = import_data("data/VolcanoPlots/" + filename)

        selected_genes = gene_hyperlinks(selected_data, ret_genes)

        figure = draw_volcano(data_dict["data"], fold_changes, selected_data, selected_genes)
        return ret_string, ret_genes, figure
                    
        return selected_genes

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

if ".DS_Store" in files:
    files.remove(".DS_Store")

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
        ]),
        html.Div([
        dcc.Markdown("""
            **Clicked Genes**
                Genes that were clicked in the volcano plot
            """),
            html.Pre(id='click-data', style=styles['pre']),
        ]),
    ]),
])
