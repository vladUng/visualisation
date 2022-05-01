from plotlyflask.plotlydash.main import app as dash_app
from plotlyflask.plotlydash.main import menu

from dash import dcc
from dash import html
from dash.dependencies import Input, Output, State
import dash_bio as dashbio
import plotly.express as px


import numpy as np
import pandas as pd

from os import path, walk
import time

pd.options.plotting.backend = "plotly"

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

### Volcano
def draw_volcano(df, fold_changes, selected_data, selected_genes):
    """ Drawing the volcano plot

    # There is something weird going on with indexes because we're setting customData to df['genes'].
    # It can get confusing as the data in df['genes'] will be put in alphabetical order and not correlated with the [x,y] points
    #  Thus the real value of the genes will be different and is passed through selected_genes

    Args:
        df (dataframe): _description_
        fold_changes (dataframe): _description_
        selected_data (list): _description_
        selected_genes (list): the list of selected genes

    Returns:
        fig: figure object from plotly
    """
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
        highlight=True, 
        annotation="group"
    )
    
    x_col = "fold_change"
    y_col = "-log10(q)"

    # show the selected points
    if selected_data:
        ranges = selected_data['range']
        selection_bounds = {'x0': ranges['x'][0], 'x1': ranges['x'][1],
                            'y0': ranges['y'][0], 'y1': ranges['y'][1]}

        # finding the selected point
        selected_points = [gene["customdata"] for gene in selected_data["points"]]
        selected_idxs = df[df["genes"].isin(selected_points)].index
        fig.update_traces(selectedpoints=selected_idxs, mode='markers', unselected={'marker': { 'opacity': 0.3 } })
    else:
        selection_bounds = {'x0': np.min(df[x_col] - 2), 'x1': np.max(df[x_col] + 2),
                            'y0': np.min(df[y_col] - 2), 'y1': np.max(df[y_col]) + 2}

    fig = show_selected_genes_vulcano(df, selected_genes, fig)

    fig.update_traces(customdata=df["genes"], textposition="bottom left")
    fig.update_layout(dragmode='select')

    fig.add_shape(dict({'type': 'rect',
                        'line': { 'width': 1, 'dash': 'dot', 'color': 'darkgrey' } },
                       **selection_bounds))
    

    return fig

def show_selected_genes_vulcano(df, selected_genes, fig):
    custom_traces = create_custom_traces(selected_genes=selected_genes)
    colors =  px.colors.qualitative.Vivid
    for idx, trace in enumerate(custom_traces): 
        fig.add_trace(create_gene_trace(df, trace["genes"], name=trace["title"], marker_color=colors[idx]))

    return fig

### Pi Plot
def draw_pi_plot(df_1, df_2, file_1, file_2, selected_data, selected_genes):

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

    x_col = "x"
    y_col = "y"

    # Add the selected points
    if selected_data:
        ranges = selected_data['range']
        selection_bounds = {'x0': ranges['x'][0], 'x1': ranges['x'][1],
                            'y0': ranges['y'][0], 'y1': ranges['y'][1]}

        # There is something weird going on with indexes
        selected_genes = [gene["customdata"][0] if "customdata" in gene.keys() else gene["text"] for gene in selected_data["points"]]

        # [gene for gene in selected_data["points"] if "customdata" not in gene.keys() else geme["text"] ]

        selected_idxs = dummy_df[dummy_df["genes"].isin(selected_genes)].index

        fig.update_traces(selectedpoints=selected_idxs,
            mode='markers', 
            unselected={'marker': { 'opacity': 0.3 }
            })
    else:
        offset = 10
        selection_bounds = {'x0': np.min(dummy_df[x_col] - offset), 'x1': np.max(dummy_df[x_col] + offset),
                            'y0': np.min(dummy_df[y_col] - offset), 'y1': np.max(dummy_df[y_col]) + offset}

    fig.update_layout(dragmode='select')
    fig.add_shape(dict({'type': 'rect',
                        'line': { 'width': 1, 'dash': 'dot', 'color': 'darkgrey' } },
                       **selection_bounds))

    fig = add_anottations(first_df, second_df, dummy_df, fig)
    fig = show_selected_genes_pi(first_df.reset_index(), second_df.reset_index(), fig, selected_genes)

    return fig

def show_selected_genes_pi(df_1, df_2, fig, selected_genes):
    custom_traces = create_custom_traces(selected_genes)
    colors =  px.colors.qualitative.Vivid
    for idx, trace in enumerate(custom_traces): 
        fig.add_trace(create_gene_trace(df_1, trace["genes"], name=trace["title"], marker_color=colors[idx], df_2=df_2))

    return fig 

def add_anottations(first_df, second_df, dummy_df, fig):
    offset = 10
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

    fig.add_annotation(showarrow=True,
                   arrowhead=1,
                   align = 'right',
                   x=first_df["x"].max() + offset, y=0,
                   text=first_df.loc[first_df["x"] == first_df["x"].max()]["comp_1"].values[0],
                   opacity=0.7)

    fig.add_annotation(showarrow=True,
                arrowhead=1,
                align = 'left',
                x=first_df["x"].min() - offset, y= 0,
                text=first_df.loc[first_df["x"] == first_df["x"].min()]["comp_1"].values[0],
                opacity=0.7)
        
    fig.add_annotation(showarrow=True,
                   arrowhead=1,
                   align = 'right',
                   y=second_df["y"].max() + offset, x=0,
                   text=second_df.loc[second_df["y"] == second_df["y"].max()]["comp_2"].values[0],
                   opacity=0.7)

    fig.add_annotation(showarrow=True,
                arrowhead=1,
                align = 'right',
                y=second_df["y"].min() - offset, x =0,
                text=second_df.loc[second_df["y"] == second_df["y"].min()]["comp_2"].values[0],
                opacity=0.7)

    return fig

### Scatter Plot

### Common plotting
def create_custom_traces(selected_genes = None):

    ifnq_genes = ['B2M', 'BATF2', 'BTN3A3', 'CD74', 'CIITA', 'CXCL10', 'CXCL9',
       'EPSTI1', 'GBP1', 'GBP4', 'GBP5', 'HAPLN3', 'HLA-DPA1', 'IDO1',
       'IFI30', 'IFI35', 'IFIT3', 'IFITM1', 'IL32', 'IRF1', 'NLRC5',
       'PARP9', 'PSMB8-AS1', 'PSMB9', 'SAMHD1', 'SECTM1', 'STAT1', 'TAP1',
       'TNFSF13B', 'TYMP', 'UBE2L6', 'WARS1', 'ZBP1']

    luminal_markers = ["KRT20", "PPARG", "FOXA1", "GATA3", "SNX31", "UPK1A", "UPK2", "FGFR3"]

    basal_markers = ["CD44", "KRT6A", "KRT5", "KRT14", "COL17A1"]

    squamos_markers = ["DSC3", "GSDMC", "TCGM1", "PI3", "TP63"]

    immune_markers = ["CD274", "PDCD1LG2", "IDO1", "CXCL11", "L1CAM", "SAA1"]

    neural_diff = ["MSI1", "PLEKHG4B", "GNG4", "PEG10", "RND2", "APLP1", "SOX2", "TUBB2B"]

    ryan_genes = ["FGFR3", "EGFR", "TP53"]
    rarg_genes = ['ADGRF1', 'ALDH3B1', 'ANKRD13D', 'ANXA11', 'B3GNT3', 'BHLHE41',
       'CAPN1', 'CDK2AP2', 'CIB1', 'CPTP', 'CTSH', 'GALE', 'HSH2D',
       'MMEL1', 'OAS1', 'PLCD3', 'PNPLA2', 'PRKCD', 'RARG', 'RNPEPL1',
       'TCIRG1', 'THEM6', 'TNFAIP2', 'UBA7', 'UNC93B1', 'VAMP8', 'VSIG2']

    # both ba_sq_inf and mes_like
    lund_qtc1 = ["FLI1", "FOXP3", "ILKZF1", "IRF4", "IRF8", "RUNX3", "SCML4", "SPI1", "STAT4", "TBX21", "TFEC"]
    lund_qtc2 = ["AEBP1", "BNC2", "GLI2", "GLIS1", "HIC1", "MSC", "PPRX1", "PPRX2", "TGFB1I1", "TWIST1"]
    lund_qtc3 = ["EBF1", "HEYL", "LEF1", "MEF2C", "TCF4", "ZEB1", "ZEB2"]
    lund_qtc8 = ["GATA5", "HAND1", "HAND2", "KLF16"]
    lund_qtc17 = ["ARID5A", "BATF3", "VENTX"]

    lund_ba_mes = lund_qtc1 + lund_qtc2 + lund_qtc3 + lund_qtc8 + lund_qtc17

    lund_ba_sq = ["BRIP1", "E2F7", "FOXM1", "ZNF367", "IRF1", "SP110", "STAT1"]
    lund_mes = ["TP53", "RB1", "FGFR3", "ANKHD1", "VIM", "ZEB2"]
    ba_sq_inf = ["CDH3", "EGFR"]
    

    custom_traces = []
    if selected_genes:
        custom_traces.append({"genes":selected_genes, "title": "Selected genes"})

    # SB
    custom_traces.append({"genes":ifnq_genes, "title": "SB_IFNQ"})
    # custom_traces.append({"genes":rarg_genes, "title": "SB_RARG"})
    # TCGA
    custom_traces.append({"genes":luminal_markers, "title": "TCGA_luminal"})
    custom_traces.append({"genes":basal_markers, "title": "TCGA_basal"})
    custom_traces.append({"genes":squamos_markers, "title": "TCGA_squamos"})
    custom_traces.append({"genes":immune_markers, "title": "TCGA_immune"})
    custom_traces.append({"genes":neural_diff, "title": "TCGA_neuroendocrine"})

    # lund
    custom_traces.append({"genes":lund_ba_mes, "title": "Mes-like + Ba/sq-inf"})
    custom_traces.append({"genes":lund_ba_sq, "title": "Lund ba/sq"})
    custom_traces.append({"genes":lund_mes, "title": "Mes-like"})
    custom_traces.append({"genes":ba_sq_inf, "title": "Ba/sq-inf"})

    return custom_traces

def create_gene_trace(df, genes, name="custom genes", marker_color="yellow", marker_size=6, df_2=None): 

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

    trace = dict(type='scatter', x=x, y= y,  showlegend=True, marker=markers, text=selected_df["genes"], mode="markers" , name=name)

    return trace

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

def create_urls(selected_data):
    """ For a given list of genes strings generate the gene cards urls

    Args:
        genes ([String]): List of genes strings

    Returns:
        [String]: URls of gene cards
    """
    base_url = "https://www.genecards.org/cgi-bin/carddisp.pl?gene="
    urls, selected_genes = [], []
    if selected_data is not None:
       
        for point in selected_data["points"]:
            gene = ""
            if "text" in point.keys():
                if "<br>" in point["text"]:
                    gene = point["text"].split("<br>")[1].split(" ")[1] 
                else:
                    gene = point["text"]
            else:
                # for pi plot
                gene = point["customdata"][0]

            selected_genes.append(gene)

        selected_genes = sorted(list(set(selected_genes)))
        urls.append(
            html.Div([
                 dcc.Textarea(
                     id="textarea_id",
                     value="[\"" + '","'.join(selected_genes) + "\"]", style={"height": 100, "width": 300},
              ),
            ])
        )
            
        # html.Div("Copy/paste list of genes: [" + '","'.join(selected_genes) + "]", style={"column-count": "2"}))
        for gene in selected_genes:
            urls.append(html.A(gene +", ", href=base_url+gene, target="_blank"))
            urls.append(html.Br())

    return urls, selected_genes

# Callbacks
def init_callbacks(dash_app): 
    @dash_app.callback(
        [Output("volcano-text-output", "children"), Output('click-data', 'children'),
         Output('figure-volcano', "figure"), Output("figure-scatter", "figure")],
        [Input('plot-volcano', 'n_clicks'),
         Input('default-volcanoplot-input', 'value')], 
          Input('figure-volcano', 'selectedData'),
        [State("select-file-volcano", "value")]
    )
    def plotVolcano(btn, fold_changes, selected_data, filename):
        ret_string = ""
        genes_div, selected_genes = [], []
        data_dict = import_data(DATA_PATH + filename)
 
        genes_div, selected_genes = create_urls(selected_data) 

        figure = draw_volcano(data_dict["data"], fold_changes, selected_data, selected_genes=selected_genes)
        return ret_string, genes_div, figure, figure
                    
    @dash_app.callback(
    [Output("pi-plot-text-output", "children"),
        Output('figure-pi-plot', "figure"), Output('click-data-pi', 'children')],
    [Input('plot-pi', 'n_clicks'),  Input('figure-pi-plot', 'selectedData')],
    [State("select-file-pi-1", "value"), State("select-file-pi-2", "value")]
    )
    def plotPi(btn, selected_data, file_1, file_2):
        ret_string = ""
        data_dict_1 = import_data(DATA_PATH + file_1)
        data_dict_2 = import_data(DATA_PATH + file_2)

        genes_div, selected_genes = create_urls(selected_data) 

        figure = draw_pi_plot(data_dict_1["data"], data_dict_2["data"], file_1, file_2, selected_data, selected_genes=selected_genes)
        return ret_string, figure, genes_div


DATA_PATH = "data/VolcanoPlots/"
files = next(walk(DATA_PATH), (None, None, []))[2]

files.sort()

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
        html.Div(id="figures", children=[
            html.Div(id="pi-vulcano-fig",  style={"column-count": "2"}, children=[
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
            html.Div(id="scatter-plot", children=[
                dcc.Graph(
                    id='figure-scatter',
                    figure={"data": {},
                            "layout": {"title": "No Scatter plot displayed", "height": 800}
                    }),                
            ])
        ]),
        html.Div(style={"column-count": "1"}, children = [
            dcc.Markdown("""
                **Pi Plot selected Genes**
                    Genes that were clicked in the pi plot
                """),
            html.Pre(id='click-data-pi', style=styles['pre']),
            dcc.Markdown("""
                **Volcano selected Genes**
                    Genes that were clicked in the volcano plot
                """),
            html.Pre(id='click-data', style=styles['pre']),
        ]),
        html.Div(style={"height":"200px"})
    ]),
])
