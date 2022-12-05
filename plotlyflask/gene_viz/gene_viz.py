import time
from datetime import datetime
from os import name, path

import dash
import numpy as np
import pandas as pd
import plotly.graph_objects as go
from dash import dash_table, dcc, html
from dash.dependencies import Input, Output, State

from plotlyflask.gene_viz import compare_genes, datasets_plot, processing
from plotlyflask.plotlydash.main import app as dash_app
from plotlyflask.plotlydash.main import menu


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
        ret_dict["all_tsv"] = pd.read_csv(
            base_path + "all_data.tsv", delimiter="\t")
        raw_df = pd.read_csv(base_path + "metadata_v3.tsv", delimiter="\t")
        # we need to do some pre-processing to duplicate the data _incl for fields
        ret_dict["metadata"], ret_dict["all_tsv"] = processing.process_data(
            raw_df, ret_dict["all_tsv"])
        print("Finished loading the data in {}".format(time.time()-start_time))
        return ret_dict
    else:
        return None

##### HTML functions #####

def create_datasets_tick_box(data_dict):
    # itereate through each of the dictionary and create the html code
    options = []
    datasets = data_dict["metadata"]["Dataset"].unique()
    for dataset in datasets:
        options.append({
            "label": dataset,
            "value": dataset,
        })
    return html.Div([
        # html.H5("Choose between dataset and subset"),
        # dcc.RadioItems( id="datasource-type",
        #     options=[
        #         {'label': 'Dataset', 'value': 'dataset'},
        #         {'label': 'Subset name', 'value': 'subset'},
        #     ],
        #     value='dataset'
        # ),
        html.Label([html.H5('Select the data sources'),
                    dcc.Link('Dataset information', href="https://docs.google.com/document/d/1yJYfFn6kdUS2KAp1tuGDLEzuR2hsscdESQTpLCxH1wU/edit")]),
        html.Br(),
        dcc.Checklist(
            id="data-checklist",
            options=options, value=datasets,
            style={"column-count": "2"}),
    ])


def create_dataset_panel(data_dict):
    # column names should be the same as the ones from the dataframe containing the gen
    return html.Div(id="gene-controls", children=[
        create_datasets_tick_box(data_dict),
        html.Br(),
        dcc.RadioItems(
            id="select-all-none",
            options=[
                {'label': 'All datasets', 'value': 'all'},
                {'label': 'Include incl_ datasets', 'value': 'incl'},
                {'label': 'Exclude incl_ datasets', 'value': 'excl'},
                {'label': 'None datasets', 'value': 'none'},
            ]
        ),
        html.Br(),
        html.Button(id='plot-gene', n_clicks=0, children='Plot'),
        html.Div(id='search-result'),
        html.Br(),
    ])


def create_gene_search():
    return html.Div(id="gene-search", children=[
        html.Hr(),
        html.H6("Enter the gene you want to analyse"),
        html.Div(["Gene name: ",
                  dcc.Input(id='gene-input', value="", type='text'),
                  html.Button(id='gene-btn', n_clicks=0, children='Add gene'),
                  html.Hr(),
                  html.Div(
                      [html.Label(["Selected gene(s)",
                                   dcc.Dropdown(
                                       id='genes-dropdown',
                                       options=[], value=[], multi=True)])
                       ],
                      style={'width': '20%', 'display': 'inline-block'}),
                  html.Div(id='display-selected-values')
                  ]),
        html.Hr(),
    ])


def metadata_menu():
    # we need to remove some non-numerical values
    x_axis_options = ["NHU_differentiation", "Tissue", "Gender",
                      "Diagnosis", "Substrate", "Dataset", "subset_name", "Sample"]
    return html.Div(id="metadata-menu", style={"column-count": "1", "margin-top": "40%"},
                    children=[
        html.H5("Figure Configuration"),
        html.H6('Indicate tight barrier? (>500)'),
        dcc.Checklist(id="ter-checklist",
                      options=[
                          {'label': 'Tight TER barrier', 'value': 'tight'}
                      ]),
        html.H6('Select the figure plot type'),
        dcc.RadioItems(
            id="select-plot-type", value='swarm',
            options=[
                {'label': 'Swarm', 'value': 'swarm', },
                {'label': 'Violin', 'value': 'violin'},
                {'label': 'Violin and Points', 'value': 'violin_points'},
                {'label': 'Box', 'value': 'box'},
                {'label': 'Box and Points', 'value': 'box_points'}
            ]
        ),
        html.H6('Select the X Axis type'),
        dcc.RadioItems(
            id="select-xaxis-type", value='linear',
            options=[
                {'label': 'Linear', 'value': 'linear'},
                {'label': 'Log10 ', 'value': 'log10'}
            ]
        ),
        html.H6('Select the X Axis'),
        dcc.Dropdown(
            id="xaxis-dropdown", value="subset_name",
            options=[{"label": i, "value": i} for i in x_axis_options]
        )
    ])


def export_panel():
    return html.Div(id="export-panel", children=[
        html.H5("Change the configuration of the plot to be saved"),
        html.Div(["Width: ",
                  dcc.Input(id='plot-width', value=1280, type='number')]),
        html.Div(["Height: ",
                  dcc.Input(id='plot-height', value=720, type='number')]),
        html.Div(["Scale: ",
                  dcc.Input(id='plot-scale', value=2, type='number')]),
        html.Br(),
        html.Div(style={"column-count": "2"}, children=[
            html.Button(id='export-plot', n_clicks=0, children='Export Plot'),
            html.Button(id='export-data', n_clicks=0, children='Export CSV')
        ]),
        html.Div(id="output-export-csv")
    ])


def init_callbacks(dash_app, data_dict):
    # latest fig to save
    @dash_app.callback(
        [Output('search-result', "children"), Output('dataset-plots-viz', "figure"), Output('intermediate-value', 'children'),
         Output('pearson-table', 'columns'), Output('pearson-table', 'data')],
        Input('plot-gene', 'n_clicks'),
        [State("genes-dropdown", "value"),  State("data-checklist", "value"), State("ter-checklist", "value"),
         State("xaxis-dropdown", "value"), State("select-plot-type", "value"), State("select-xaxis-type", "value")]
    )
    def button_router(btn_1, user_input, datasets_selected, ter, xaxis, plot_type, xaxis_type):

        ret_data = pd.DataFrame().to_json(date_format='iso', orient='split')
        if len(user_input):
            all_tsv, metadata = data_dict["all_tsv"], data_dict["metadata"]
            if len(user_input) == 1:
                gene = user_input[0].strip()
                # search for the gene
                found_in = all_tsv[all_tsv["genes"].str.fullmatch( "\\b"+gene+"\\b", case=False, na=False)]
                # create the figure
                ret_string, last_fig, prcsd_data = datasets_plot.plot(
                    found_in, metadata, gene, datasets_selected, ter, plot_type, xaxis, xaxis_type)
                # transformed the filtered data to json
                ret_data = prcsd_data.to_json(
                    date_format='iso', orient='split')

                return ret_string, last_fig, ret_data, [], []
            else:
                # we have multiples genes to look for
                first_gene = all_tsv[all_tsv["genes"].str.fullmatch("\\b"+ user_input[0].strip() +"\\b", case=False, na=False)]
                second_gene = all_tsv[all_tsv["genes"].str.fullmatch("\\b"+ user_input[1].strip( )+"\\b", case=False, na=False)]
                found_in = pd.concat([first_gene, second_gene])

                # create the figure
                ret_string, last_fig, prcsd_data = compare_genes.plot( found_in, metadata, user_input, datasets_selected, plot_type, xaxis_type)

                # normal
                gene_A, gene_B = found_in["genes"].values[0], found_in["genes"].values[1]
                dummy_df = found_in.drop("genes", axis=1).transpose()
                dummy_df.columns = found_in["genes"].values
                pearson_corr, spearman_corr = processing.correlations(dummy_df, gene_A, gene_B)

                # for log10
                dummy_df = pd.DataFrame()
                dummy_df[found_in["genes"].values[0]] = np.log10(list(found_in.iloc[0, 1:].values + 1))
                dummy_df[found_in["genes"].values[1]] = np.log10(list(found_in.iloc[1, 1:].values + 1))
                # we don't need to calculate for spearman, but it's easier to the function
                pearson_corr_log, _ = processing.correlations( dummy_df, gene_A, gene_B, isLog=True)

                corelation_cols = [{"id": "corr_metric", "name": "Correlation Coefficient"},
                                   {"id": "value", "name": "{} vs {}".format(gene_A, gene_B)}]

                table_data = [
                    {"corr_metric": "Pearson Correlation",
                        "value": pearson_corr[0]["gene_b"]},
                    {"corr_metric": "Pearson Correlation log10(TPM+1)",
                     "value": pearson_corr_log[0]["gene_b"]},
                    {"corr_metric": "Spearman Correlation",
                        "value": spearman_corr[0]["gene_b"]}
                ]

                ret_data = prcsd_data.to_json(date_format='iso', orient='split')
                return ret_string, last_fig, ret_data, corelation_cols, table_data

        return "Search for gene and click submit", {"data": {}}, ret_data, [], []

    @dash_app.callback(
        [Output('genes-dropdown', 'options'), Output('genes-dropdown', 'value'),
         Output('display-selected-values', 'children')],
        [Input("gene-btn", 'n_clicks'), Input("gene-input", "n_submit")],
        [State("gene-input", "value"), State("genes-dropdown",
                                             "options"), State("genes-dropdown", "value")]
    )
    def update_gene_dropdown(btn, n_submit, new_gene, genes_options, selected_genes):
        retMessage = ""
        if new_gene != '':
            all_tsv, metadata = data_dict["all_tsv"], data_dict["metadata"]
            found_in = all_tsv[all_tsv["genes"].str.fullmatch(
                "\\b"+new_gene+"\\b", case=False, na=False)]

            if not found_in.empty:
                # check if duplicated
                geneExists = False
                for option in genes_options:
                    if new_gene in option.values():
                        geneExists = True
                        break

                if not geneExists:
                    genes_options.append(
                        {"label": new_gene, "value": new_gene})
                    retMessage = "Added gene {}".format(new_gene)

                    if len(selected_genes) > 1:
                        retMessage = "You can only select 2 genes, please remove one."
                    else:
                        selected_genes.append(new_gene)
                else:
                    if len(selected_genes) < 2 and not geneExists:
                        selected_genes.append(new_gene)
                    retMessage = "Gene {} is already in the options".format(
                        new_gene)
            else:
                retMessage = "We couldn't find the gene {} in the selected dataset".format(
                    new_gene)

        return genes_options, selected_genes, retMessage

    @dash_app.callback(
        Output("data-checklist", "value"),
        [Input("select-all-none", "value")],
        [State("data-checklist", "value")]
    )
    def select_all_none_dasets(radio_value, datasets):
        if radio_value == "none":
            return []
        elif radio_value == "excl":
            incl_cols = [
                col for col in data_dict["metadata"].columns if "incl_" in col]
            for incl_col in incl_cols:
                if incl_col in datasets:
                    datasets.remove(incl_col)
            return datasets
        elif radio_value == "incl":
            incl_cols = [
                col for col in data_dict["metadata"].columns if "incl_" in col]
            # include only the incl_cols that were not selected
            datasets += list(set(incl_cols) - set(datasets))
            return datasets
        else:
            return list(data_dict["metadata"]["Dataset"].unique())

    @dash_app.callback(
        Output("output-export-csv", 'children'),
        [Input("export-plot", "n_clicks"), Input("export-data",
                                                 "n_clicks"), Input('intermediate-value', 'children')],
        [State("genes-dropdown", "value"), State("plot-width", "value"), State("plot-height", "value"),
         State("plot-scale", "value"), State("dataset-plots-viz", "figure"), State("data-checklist", "value")]
    )
    def export_plot(btn_1, btn_2,  data_json, user_input, width, height, scale, figure, datasets):
        # checking os
        if name == 'nt':
            path = "exported_data\\"
        else:
            path = "./exported_data/"
        ret_string = "There is no figure to save"

        ctx = dash.callback_context
        # get the button_id
        if not ctx.triggered:
            button_id = 'No clicks yet'
        else:
            button_id = ctx.triggered[0]['prop_id'].split('.')[0]

        if button_id == "export-plot":
            if figure["data"]:
                fig = go.Figure(figure)
                filepath = "{}plot-jbu-viz-{}{}".format(
                    path, datetime.now().strftime('%d-%m-%Y--%H-%M-%S'), ".pdf")
                fig.write_image(filepath, width=width,
                                height=height, scale=scale)
                ret_string = "Plot saved to {}".format(filepath)
        elif button_id == "export-data":
            if len(user_input) == 1:
                figure_data_df = pd.read_json(data_json, orient='split')
                filepath = "{}data-jbu-viz-{}{}".format(
                    path, datetime.now().strftime('%d-%m-%Y--%H-%M-%S'), ".csv")
                filepath_all = "{}all-data-jbu-viz-{}{}".format(
                    path, datetime.now().strftime('%d-%m-%Y--%H-%M-%S'), ".csv")
                # we need to remove the duplicated samples
                unique_subsets = figure_data_df["subset_name"].unique()
                export_df = pd.DataFrame()
                for subset_name in unique_subsets:
                    dummy_df = figure_data_df[figure_data_df["subset_name"]
                                              == subset_name]["TPM"].T
                    # we need to rename the series so that has a unique name, otherwise the below line won't work
                    export_df = pd.concat(
                        [export_df, dummy_df.rename(subset_name)], axis=1)

                # we need to put the nan values at the bottom and to that we need to sort each column individualy and re-concat
                export_df = pd.concat([export_df[col].sort_values().reset_index(
                    drop=True) for col in export_df], axis=1, ignore_index=True)
                export_df.columns = [
                    subset+"_TPM" for subset in unique_subsets]

                export_df.to_csv(filepath, index=False)
                all_data = figure_data_df[figure_data_df["subset_name"].isin(
                    unique_subsets)]
                all_data = all_data.dropna(axis=1, how='all')

                # remove _incl from samples and datasets
                samples = []
                for _, row in all_data.iterrows():
                    sample = row["Sample"]
                    if "incl_" in sample:
                        samples.append("-".join(sample.split("_")[:-2]))
                    else:
                        samples.append(sample)

                all_data["Sample"] = samples

                tpm_col = all_data.pop("TPM")
                all_data.insert(1, "TPM", tpm_col)
                all_data.to_csv(filepath_all, index=False)

                ret_string = "Data saved to {} and with the metadata to: {}".format(
                    filepath, filepath_all)

            else:
                figure_data_df = pd.read_json(data_json, orient='split')
                filepath = "{}data-jbu-viz-{}{}".format(
                    path, datetime.now().strftime('%d-%m-%Y--%H-%M-%S'), ".csv")

                # sample, subset_name and the two genes (last ones)
                selected_cols = np.append(
                    figure_data_df.columns[-2:].values, figure_data_df.columns[:2])
                # we need to remove the duplicated samples
                figure_data_df[selected_cols].to_csv(filepath, index=False)
                ret_string = "Data saved to {}".format(filepath)

        return ret_string

    """Create a Plotly Dash dashboard."""


# update_all_datasets("data/")
data_dict = import_data("JBU_data/")
init_callbacks(dash_app, data_dict)

# Create Layout
layout = html.Div(children=[
    menu,
    html.H1(children='JBU visualisation tool'),
    html.Div(id="parrent-all-viz", children=[
        create_gene_search(),
        html.Div(id="parent-controllers-viz", style={"column-count": "2", "width": "fit-content"},
                    children=[
            create_dataset_panel(data_dict),
            metadata_menu(),
        ]),
        dcc.Graph(
            id='dataset-plots-viz',
            figure={"data": {},
                    "layout": {"title": "No Gene Search", "height": 300}
                    }),
        html.Div([
            html.H5("Correlations Table"),
            html.Hr(),
            dash_table.DataTable(id="pearson-table", columns=[], data=[],
                      style_table={"width": "30%"},
                      style_header={
                'backgroundColor': 'rgb(230, 230, 230)', 'fontWeight': 'bold'},
                style_cell={'textAlign': 'left'}),
            html.Hr(),
            html.Div(["Note: Pearson (raw and log10 transformed) and Spearman correlation values computed here are for guidance only. Values will be inaccurate/inappropriate with major outliers. Pearson correlations will be incorrect when data is non-normally distributed - common with expression data, hence the log10(TPM+1) transformation. Spearman is a non-parametric test built on rank order, so will not be affected by log transformations. Data should be downloaded (button below) and inspected before presentation/interpretation."]),
            html.Hr()]),
        export_panel(),
        html.Div(id='intermediate-value', style={'display': 'none'})
    ]),
])
