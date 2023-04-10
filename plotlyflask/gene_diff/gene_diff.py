#!/usr/bin/env python
# -*-coding:utf-8 -*-
'''
@File    :   gene_diff.py
@Time    :   2022/11/11 13:11:13
@Author  :   Vlad Ungureanu
@Version :   1.0
@Contact :   vlad.ungureanu@york.ac.uk
@Desc    :   Main file for Differential Expression Analysis page on the JBU visualisation tool.

The purpose of this fail is to create the layout for the page (see the global variable layout) and call the plot functions. Currently the following features are supported:

1. Volcano plot (volcano.py)
2. Pi plot (pi.py)
3. Scatter plot (scatter.py)

At the end of the pages the Gene cards and arrays of the selected genes are shown.

'''


from os import walk

from dash import dcc, html
from dash.dependencies import Input, Output, State

# Import the plots
from plotlyflask.gene_diff import pi, scatter, shared, volcano
from plotlyflask.plotlydash.main import app as dash_app
from plotlyflask.plotlydash.main import menu

styles = {
    'pre': {
        'border': 'thin lightgrey solid',
        'overflowX': 'scroll'
    }
}


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
            data_dict = shared.import_data(DATA_PATH + filename)
            prev_filename = filename
 
        genes_div, selected_genes = create_urls(selected_data) 

        vulcano_fig = volcano.draw_volcano(data_dict["data"], fold_changes, selected_data, selected_genes=selected_genes)
        scatter_fig = scatter.draw_scatter(data_dict["data"], selected_genes)

        return ret_string, genes_div, vulcano_fig, scatter_fig

    @dash_app.callback(
        Output('click-data-scatter', 'children'),
        Input('figure-scatter', 'selectedData'),
    )
    def scatterSelection(selected_data):

        genes_div, _ = create_urls(selected_data) 
        return genes_div

    @dash_app.callback(
    [Output("pi-plot-text-output", "children"),
        Output('figure-pi-plot', "figure"), Output('click-data-pi', 'children')],
    [Input('plot-pi', 'n_clicks'),  Input('figure-pi-plot', 'selectedData')],
    [State("select-file-pi-1", "value"), State("select-file-pi-2", "value")]
    )
    def plotPi(btn, selected_data, file_1, file_2):
        ret_string = ""
        data_dict_1 = shared.import_data(DATA_PATH + file_1)
        data_dict_2 = shared.import_data(DATA_PATH + file_2)

        genes_div, selected_genes = create_urls(selected_data) 

        figure = pi.draw_pi_plot(data_dict_1["data"], data_dict_2["data"], file_1, file_2, selected_data, selected_genes=selected_genes)
        return ret_string, figure, genes_div

# Get all the files
DATA_PATH = "JBU_data/VolcanoPlots/"
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
    html.Div(style={"column-count": "2"}, children = [
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