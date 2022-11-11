from plotlyflask.plotlydash.main import app as dash_app
from plotlyflask.plotlydash.main import menu

from dash import dcc
from dash import html
from dash.dependencies import Input, Output, State
import dash_bio as dashbio
import plotly.express as px
import plotly.graph_objects as go

import numpy as np
import pandas as pd

from os import path, walk
import time

from plotlyflask.utilities import scatter_plot as sp
from plotlyflask.utilities import gene_markers as gm

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
        point_size=10,
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
        fig.update_traces(marker=dict(size=10, opacity=0.4), selector=dict(mode='markers'))


    fig = show_selected_genes_vulcano(df, selected_genes, fig)

    fig.update_traces(customdata=df["genes"], textposition="top right")
    fig.update_layout(dragmode='select')

    fig.add_shape(dict({'type': 'rect',
                        'line': { 'width': 1, 'dash': 'dot', 'color': 'darkgrey' } },
                       **selection_bounds))

    return fig

def show_selected_genes_vulcano(df, selected_genes, fig):
    custom_traces = create_custom_traces(selected_genes=selected_genes)
    colors =  px.colors.qualitative.Bold + px.colors.qualitative.Vivid 
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
    
    dummy_df["main_colour"] = "PI_plot"
    # show in the hover everything apart from the main_colour column
    fig = px.scatter(dummy_df, x="x", y="y", hover_data=dummy_df.columns, color="main_colour", title=title)

    fig.update_traces(marker=dict(size=10, opacity=0.4), selector=dict(mode='markers'))
    
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
            mode='markers+text', 
            unselected={'marker': { 'opacity': 0.3 }
            })
    
        fig.add_shape(dict({'type': 'rect',
                        'line': { 'width': 1, 'dash': 'dot', 'color': 'darkgrey' } },
                       **selection_bounds))
    else:
        offset = 10
        selection_bounds = {'x0': np.min(dummy_df[x_col] - offset), 'x1': np.max(dummy_df[x_col] + offset),'y0': np.min(dummy_df[y_col] - offset), 'y1': np.max(dummy_df[y_col]) + offset}

    fig = add_anottations(first_df, second_df, dummy_df, fig)
    fig = show_selected_genes_pi(first_df.reset_index(), second_df.reset_index(), fig, selected_genes)

    return fig

def show_selected_genes_pi(df_1, df_2, fig, selected_genes):
    custom_traces = create_custom_traces(selected_genes)
    colors =  px.colors.qualitative.Bold + px.colors.qualitative.Vivid 
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
def plt_fc_lines(fig, df, fc_values):

    # max_value = round(df[df.columns.values[1:-1]].max().max())
    # get the maximum values across fold_change and the cluster groups
    med_cols = [col for col in df.columns if "med" in col and "fold" not in col]
    groups = ["fold_change"] + list(med_cols)
    max_value = round(df[groups].max().max())

    colours = ["Black", "Green", "goldenrod", "Red", "orange", "Green", "Purple"]
    line_type = [None, None, "dot", "dash", "dash"]
    for value in fc_values:

        if value:
            color = colours[value]
        else:
            color = "Black"

        line_dict = dict(color=color, width=3,  dash=line_type[value])

        fig.add_shape(type='line',
                x0=value, y0=0,
                x1=max_value, y1=max_value-value,
                line=line_dict,
                xref='x',yref='y', name="log2(FC)={}".format(value))

        if value:
            color = colours[value]
        else:
            color = "Black"

        line_dict = dict(color=color, width=3,  dash=line_type[value])

        fig.add_shape(type='line',
                x0=0, y0=value,
                x1=max_value-value+1, y1=max_value+1,
                line=line_dict,
                xref='x',yref='y', name="log2(FC)={}".format(value))

    return fig

def plt_scatter(df, selected_points, known_markers=False):

    # select the columns with fold change based on the median
    clusters = [col for col in df.columns if ("med" in col) and ("fold" not in col)]
    fig = px.scatter(df, x=clusters[0], y=clusters[1], hover_data=df.columns, color="significance",  color_discrete_sequence=['#636EFA', "grey", '#EF553B'])

    fig.update_traces(marker=dict(size=10, opacity=0.4), selector=dict(mode='markers'))
    
    fig = plt_fc_lines(fig, df, fc_values = [0, 1, 2, 4])

    if not selected_points.empty:
        fig.add_trace(go.Scatter(x=selected_points[clusters[0]], y=selected_points[clusters[1]], mode="markers+text", text=selected_points["genes"], hoverinfo='all', textposition="top right", name="Selected Points"))

    # # global clustered_genes
    # if not clustered_genes.empty:
    #     dmy_df = df[df["genes"].isin(clustered_genes)]
    #     fig.add_trace(go.Scatter(x=dmy_df[clusters[0]], y=dmy_df[clusters[1]],mode="markers", text=dmy_df["genes"], hoverinfo='all', textposition="top right", name="Clustered genes"))

    if known_markers:
        custom_traces = create_custom_traces()
        colors =  px.colors.qualitative.Vivid + px.colors.qualitative.Bold
        marker_size = 12
        for idx, trace in enumerate(custom_traces): 
            selected_df = df[df["genes"].isin(trace["genes"])]

            markers = {"size": marker_size, "color": selected_df.shape[0] * [colors[idx]], "symbol": "x"}
            trace = dict(type='scatter', x=selected_df[clusters[0]], y=selected_df[clusters[1]],  showlegend=True, marker=markers, text=selected_df["genes"], mode="markers+text" , name=trace["title"],  textposition="top right")

            fig.add_trace(trace)

    fig.update_layout(
        xaxis = dict(
            tickmode = 'array',
            tickvals = np.log2([1, 10, 100, 1000, 10000]),
            ticktext = ["0", "1", "10", "100", "1000", "10000"]
        ),
        yaxis = dict(
            tickmode = 'array',
            tickvals = np.log2([1, 10, 100, 1000, 10000]),
            ticktext = ["0", "1", "10", "100", "1000", "10000"]
        )
    )
    
    return fig

def draw_scatter(df, selected_genes):
    # exp = "_".join(filename.split("_")[:-2])
    # fold_change, _ = sp.prc_fc_comp(tcga_tpm_df, df, exp, mapping_cols, drop_outliers=False)

    df = sp.add_sig(df)
    selected_points = df[df["genes"].isin(selected_genes)]

    fig = plt_scatter(df, selected_points, True)
    fig.update_layout(clickmode='event+select')

    return fig

### Common plotting
def create_custom_traces(selected_genes = None):
    custom_traces = []
    if selected_genes:
        custom_traces.append({"genes":selected_genes, "title": "Selected genes"})

    ifnq_genes = ["INFG", 'B2M', 'BATF2', 'BTN3A3', 'CD74', 'CIITA', 'CXCL10', 'CXCL9','EPSTI1', 'GBP1', 'GBP4', 'GBP5', 'HAPLN3', 'HLA-DPA1', 'IDO1', 'IFI30', 'IFI35', 'IFIT3', 'IFITM1', 'IL32', 'IRF1', 'NLRC5', 'PARP9', 'PSMB8-AS1', 'PSMB9', 'SAMHD1', 'SECTM1', 'STAT1', 'TAP1','TNFSF13B', 'TYMP', 'UBE2L6', 'WARS1', 'ZBP1']

    diff_neuronal = ['ST18', 'DPYSL5', 'DCX', 'TMEM145', 'ELAVL3', 'TOP2A', 'SCN3A', 'ASXL3', 'SEZ6', 'ATP1A3', 'CELF3', 'PEX5L', 'PCSK2', 'PROX1', 'EPHA7', 'CHGB', 'UNC13A', 'RIMS2', 'CAMK2B', 'NKX2-2', 'AP3B2']

    # diff when Consensus vs Mixed and LumInf vs Mixed
    ne_dif = ["ARHGAP33","CBX5","CCNE2","CDCA5","CDCA8","CDK5R1","CDKAL1","CENPF","CLSPN","CNIH2","COCH","DCX","DPYSL5","E2F7","ENHO","FBXO43","FBXO5","FSD1","GNAZ","GPRIN1","GTSE1","HES6","HMGB2","IFT81","KCNH3","LHX2","LMNB1","MAD2L1","MCM2","MSANTD3-TMEFF1","MT3","PBK","PDXP","PIMREG","PPM1D","PRAME","PSIP1","PSRC1","PTTG1","QSOX2","RCOR2","SALL2","SCML2","SMC2","SRSF12","STMN1","TEDC2","TLCD3B","TMSB15A","TUBA1B","TUBB2B","TXNDC16","TYMS","USP1","WASF1","ZNF711", "C12orf75","CBX2","CDC25C","DEPDC1","EZH2","GAS2L3","KIF18B","NUSAP1","ODC1"]

    mixed_lumInf_diff =["ADGRF4","ADIRF","ALS2CL","ANXA8","ANXA8L1","B3GNT3","BATF","BCAS1","C10orf99","C19orf33","CH17-360D5.2","CLDN1","CTSH","CXCL17","EPS8L1","EVPL","FAM110C","FXYD3","GBP2","GNA15","GPR87","IL1RN","ITGB4","ITGB6","KCNN4","KRT19","KRT7","MPZL2","MYOF","P2RY2","PLEKHN1","PRSS22","PSCA","PTGES","PTK6","S100A11","S100A6","S100P","SDC1","SNCG","SYT8","SYTL1","TACSTD2","TINAGL1","TMEM40","UPK1A","UPK2","UPK3A","UPK3B","VGLL1","WNT7B", "ACSF2","ARHGAP27","BICDL2","CAPS","CARD11","CBLC","CLDN4","CSTB","CYP4B1","DENND2D","DTX4","EHF","ELF3","EPHA1","EPN3","FBP1","FOXQ1","GATA3","GDF15","GGT6","GPR160","GRHL1","GRHL3","IQANK1","KRT7-AS","LLNLR-231D4.1","LPAR5","METRNL","NECTIN4","OVOL1","PLA2G2F","PLA2G4F","PPL","PROM2","PSD4","RAB25","RBBP8NL","RBM47","RP4-568C11.4","S100A14","SCNN1B","SEMA4A","SPINK1","SSH3","TFAP2C","TJP3","TMC6","TMPRSS2","TRPM4","UGT1A6","VAMP8","VSIG2"]

    jens_work = ["GJB1"]

    basal_small_specific = ['AC005301.9',   'AIF1', 'APOC2', 'CALHM6', 'CCL4', 'CD163',  'CLCF1',  'DERL3',  'ELN',  'FBLN2',  'FCGR2A',  'FOLR2',  'FST',  'GZMB',  'HCST',  'IFIT3',  'ITGB2',  'KRT80',  'MFAP4',  'MMP23B',  'MS4A6A',  'MT1M',  'MZB1',  'NKG7',  'NR4A1AS',  'OLFML2A',  'OLFML3',  'RGS2',  'RP11-54H7.4',  'TRAC']

    custom_traces.append({"genes": basal_small_specific, "title": "Small Ba/Sq specific"})

    # custom_traces.append({"genes": jens_work, "title": "Jen's work"})

    # custom_traces.append({"genes":ne_dif, "title": "Diff for NE"})

    # custom_traces.append({"genes":mixed_lumInf_diff, "title": "Diff for Mixed/LumInf"})


    # SB
    # custom_traces.append({"genes":ifnq_genes, "title": "SB_IFNQ"})
    # custom_traces.append({"genes":diff_neuronal, "title": "Diff old vs remap"})

    # ryan_genes = ["FGFR3", "EGFR", "TP53"]
    custom_traces = gm.add_tcga_markers(custom_traces=custom_traces)
    # custom_traces = gm.add_lund_markers(custom_traces=custom_traces)

    # custom_traces = gm.lumInf_mixed(custom_traces=custom_traces)

    custom_traces = gm.significant_genes(custom_traces=custom_traces)

    custom_traces = gm.low_significant_genes(custom_traces=custom_traces)
    
    return custom_traces
    # return []

def create_gene_trace(df, genes, name="custom genes", marker_color="yellow", marker_size=12, df_2=None): 

    selected_df = df[df["genes"].isin(genes)]

    x, y = [], []
    if df_2 is None:
        y = -np.log10(selected_df["q"])
        x = selected_df['fold_change']
    else:
        selected_df_2 = df_2[df_2["genes"].isin(genes)]
        
        x = -np.log10(selected_df["q"]) * selected_df["fold_change"]
        y = -np.log10(selected_df_2["q"]) * selected_df_2["fold_change"]

    markers = {"size": marker_size, "color": selected_df.shape[0] * [marker_color], "symbol": "x" }

    trace = dict(type='scatter', x=x, y= y,  showlegend=True, marker=markers, text=selected_df["genes"], mode="markers+text" , name=name, textposition= 'top center')

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
        urls.append(
             html.Div([
                 dcc.Textarea(
                     id="textarea_id_2",
                     value=','.join(selected_genes), style={"height": 100, "width": 300},
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

    # optimisation for volcano
    global prev_filename, data_dict 
    prev_filename = ""

    @dash_app.callback(
        [Output("volcano-text-output", "children"), Output('click-data', 'children'),
         Output('figure-volcano', "figure"), Output("figure-scatter", "figure")],
        [Input('plot-volcano', 'n_clicks'),
         Input('default-volcanoplot-input', 'value'), 
          Input('figure-volcano', 'selectedData')],
        [State("select-file-volcano", "value")]
    )
    def plotVolcano(btn, fold_changes, selected_data, filename):
        ret_string = ""
        genes_div, selected_genes = [], []
        global prev_filename, data_dict
        if prev_filename != filename:
            data_dict = import_data(DATA_PATH + filename)
            prev_filename = filename
 
        genes_div, selected_genes = create_urls(selected_data) 

        vulcano = draw_volcano(data_dict["data"], fold_changes, selected_data, selected_genes=selected_genes)
        scatter = draw_scatter(data_dict["data"], selected_genes)

        return ret_string, genes_div, vulcano, scatter

    @dash_app.callback(
        Output('click-data-scatter', 'children'),
        Input('figure-scatter', 'selectedData'),
    )
    def scatterSelection(selected_data):

        genes_div, selected_genes = create_urls(selected_data) 
        return genes_div

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

# Get all the files
DATA_PATH = "data/VolcanoPlots/"
files = next(walk(DATA_PATH), (None, None, []))[2]

if len(files):
    files.sort()

    if ".DS_Store" in files:
        files.remove(".DS_Store")
else:
    # ToDo: Alert to browswer 
    print("There are no files available");

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
    ]),
    html.Div(id="figures", style={"column-count": "2"}, children=[
        html.Div(children=[
            html.Div(id='volcano-text-output'),
            dcc.Graph(
                id='figure-volcano',
                figure={"data": {},
                        "layout": {"title": "No Volcano displayed", "height": 800}
                        }),
        ]),
        html.Div(children=[
            html.Div(id='pi-plot-text-output'),
            dcc.Graph(
                id='figure-pi-plot',
                figure={"data": {},
                        "layout": {"title": "No π plot displayed", "height": 800}
                }),      
        ]),
    ]),
        html.Div(id="scatter-plot", children=[
            dcc.Graph(
                id='figure-scatter',
                figure={"data": {},
                        "layout": {"title": "No Scatter plot displayed", "height": 600}
                }), 
            dcc.Markdown("""
            **Scatter plot selected Genes**
                Genes that were clicked in the scatter plot plot
            """),
            html.Pre(id='click-data-scatter', style=styles['pre'])   
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
    html.Div(style={"height":"200px"}),
])



