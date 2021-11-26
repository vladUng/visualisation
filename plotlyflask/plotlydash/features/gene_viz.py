from plotlyflask.plotlydash.main import menu
from plotlyflask.plotlydash.main import app as dash_app
import numpy as np
import pandas as pd
import dash
from dash_table import DataTable
import dash_html_components as html
import dash_core_components as dcc
from dash.dependencies import Input, Output, State
import plotly.express as px
import plotly.graph_objects as go

from os import path, name
import time
from datetime import datetime

##### Functions for .tsv handling #####


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
        ret_dict["metadata"], ret_dict["all_tsv"] = process_data(
            raw_df, ret_dict["all_tsv"])
        print("Finished loading the data in {}".format(time.time()-start_time))
        return ret_dict
    else:
        return None


def update_all_datasets(base_path):
    """ 
    Function which takes the names of all the files from the master.tsv files and then updates the all tsv file

    Args:
        base_path ([String]): Path to /data folder 

    Returns:
        [type]: [description]
    """
    start_time = time.time()
    # open the master .csv file, containing the contents of the folder
    master_cv = pd.read_csv(base_path + "master.tsv", delimiter="\t")

    # iterate through the csv file and create a list of dictionaries with the following: <key> - name of the file, <value> object which contains a dataframe and description from the .csv file.
    all_tsv = pd.DataFrame()
    ret_dict = {}
    for _, row in master_cv.iterrows():
        filename, _ = row["Filename"], row["Description"]
        df = pd.read_csv(base_path+filename, delimiter="\t")
        print("Dataset: {} Shape {}".format(filename, df.shape))

        if "metadata" not in filename:
            if not all_tsv.empty:
                all_tsv = pd.merge(all_tsv, df, how="outer")
            else:
                all_tsv = df

    print("Finished loading the data in {}".format(time.time()-start_time))
    ret_dict["all_tsv"] = all_tsv

    start_time = time.time()
    all_tsv.to_csv(base_path + "all_data.tsv",  sep='\t', index=False)
    print("Update .tsv file with the new data in {}".format(
        time.time() - start_time))


def process_data(df_metadata, all_tsv):
    """ Process the metadata and the all_tsv file. The are two main operations happening:
    1. Replacing in specified cols the nan values and "?" with Not-known
    2. Fakes the _incl columns with values of "Y" to be as datasets. This involved to duplicate both the samples in metadata and in the `all_data.tsv`. We've marked the samples that are duplicated by the following rule "_incl_" + "including column initials". 
    Args:
        df_metadata ([DataFrame]): The medata
        all_tsv ([DataFrame]): TPM dataframe
    Returns:
        (DataFrame, DataFrame): The two dataframes updated
    """

    df_metadata = df_metadata.replace("?", np.nan)
    # cols to replace nan with NA string
    cols = ["NHU_differentiation", "Tissue",
            "Gender", "Diagnosis", "Substrate"]
    for col in cols:
        df_metadata[col] = df_metadata[col].fillna("NA")

    # add the incl columns to the datasets and duplicate the samples
    incl_cols = [col for col in df_metadata.columns if "incl_" in col]
    new_rows_metadata = pd.DataFrame()
    new_cols_all_tsv = pd.DataFrame()
    for incl_col in incl_cols:
        # select samples
        metadata_selected = filter_data(df_metadata, incl_col, ["Y"])
        selected_samples = metadata_selected["Sample"].values
        # get initials of the cols
        initials = "_incl_"
        for initial in incl_col.split("_"):
            initials += initial[0]
        # update the metadata
        duplicated_samples = [
            sample+initials for sample in metadata_selected["Sample"].values]
        metadata_selected = metadata_selected.drop("Sample", axis=1)
        metadata_selected["Sample"] = duplicated_samples
        metadata_selected["Dataset"] = incl_col
        metadata_selected["subset_name"] = incl_col
        new_rows_metadata = new_rows_metadata.append(
            metadata_selected, ignore_index=True)

        # update all_tsv
        all_tsv_selected = all_tsv[selected_samples]
        all_tsv_selected.columns = duplicated_samples
        new_cols_all_tsv = pd.concat(
            [new_cols_all_tsv, all_tsv_selected], axis=1)

    df_metadata = pd.merge(df_metadata, new_rows_metadata, how="outer")
    df_metadata["isTER"] = True
    for subtype in df_metadata["subset_name"].unique():
        if df_metadata[df_metadata["subset_name"] == subtype]["TER"].isnull().all():
            df_metadata.loc[df_metadata["subset_name"]
                            == subtype, "isTER"] = False

    all_tsv = pd.concat([all_tsv, new_cols_all_tsv], axis=1)
    return df_metadata, all_tsv

##### Utilities functions #####


def filter_data(df, col, values_to_keep):
    if values_to_keep:
        df = df[df[col].isin(values_to_keep)]
    return df


def sync_df(df_metadata, df):
    """ The role of this functions is to make sure that metadata df and all_data corresponds to the same samples after we've filter out some off the samples.

    Args:
        df_metadata ([DataFrame]): The medata
        df ([DataFrame]): TPM dataframe

    Returns:
        DataFrame: The metadata sorted
    """
    # the df with common samples from both df metadata and the df with tpms
    df_common = filter_data(df, "Sample", list(df_metadata["Sample"].values))

    # for cases when some samples are in a database but not others
    metadata_common = filter_data(
        df_metadata, "Sample", list(df["Sample"].values))

    # sort both by samples so that there are syncd
    metadata_sorted = metadata_common.sort_values("Sample")
    df_common = df_common.sort_values("Sample")
    metadata_sorted["TPM"] = df_common["TPM"].values
    return metadata_sorted


def samples_to_remove(df):
    # used to remove _incl samples when exporting the .csv files
    ret_samples = []
    # samples duplicated by tpm
    duplicated_samples = df[df["TPM"].duplicated()]["Sample"].values
    # but we need to check if these are the duplicates of the _incl (maybe in the future there will be samples that have the exact TPM value)
    ret_samples = [
        sample for sample in duplicated_samples if "_incl_" in sample]
    return ret_samples


def compute_correlations(df, gene_A,  gene_B, isLog=False):
    spearman_df = df.corr(method="spearman").round(6)
    pearson_df = df.corr(method="pearson").round(6)

    pearson_data, spearman_data = [], []

    for idx, row in pearson_df.iterrows():
        if isLog:
            pearson_data.append({"corr_metric": "Log10({} + 1)".format(idx),
                                 "gene_a": row[gene_A], "gene_b": row[gene_B]})
        else:
            pearson_data.append(
                {"corr_metric": idx, "gene_a": row[gene_A], "gene_b": row[gene_B]})

    for idx, row in spearman_df.iterrows():
        if isLog:
            spearman_data.append({"corr_metric": "Log10({} + 1)".format(idx),
                                 "gene_a": row[gene_A], "gene_b": row[gene_B]})
        else:
            spearman_data.append(
                {"corr_metric": idx, "gene_a": row[gene_A], "gene_b": row[gene_B]})

    return pearson_data, spearman_data

##### Plotting functions #####


def create_datasets_plot(found_in, metadata, user_input, datasets_selected, ter, plot_type, xaxis, xaxis_type):
    if not found_in.empty:
        ret_str = 'Gene {} was found'.format(user_input)

        # create the figure layout
        layout = go.Layout(
            title="Swarm plot for {}".format(found_in["genes"].values[0]),
            height=700
        )
        # prepare datafrarme
        found_in = found_in.drop(["genes"], axis=1).T.reset_index()
        found_in.columns = ["Sample", "TPM"]

        # process metadata and add to the df
        metadata_selected = filter_data(metadata, "Dataset", datasets_selected)
        processed_df = sync_df(metadata_selected, found_in)

        hover_data = ["Sample",  "Dataset", "subset_name", "shared_num_same_col",
                      "Tissue", "NHU_differentiation", "Gender", "TER", "Substrate"]
        if not processed_df.empty:
            if not datasets_selected:
                return "Select a dataset first", {"data": {}}, found_in

            config = {"x": xaxis, "y": "TPM", "color": "Dataset",
                      "hover_data": hover_data, "xaxis_type": xaxis_type, }

            # create the main trace
            fig = select_plot_type(processed_df, config,
                                   plot_type, isComparison=False, ter=ter)

            if xaxis_type == "linear":
                if processed_df["TPM"].values.max() <= 25:
                    fig.update_yaxes(range=[0, 25])
                else:
                    offset_y = processed_df["TPM"].values.max() * 0.015
                    fig.update_yaxes(
                        range=[-offset_y, processed_df["TPM"].values.max() + offset_y])

            fig.update_layout(layout)
            return ret_str, fig, processed_df

    return "Gene {} not found in any of the datasets".format(user_input), {"data": {}}, found_in


def compare_genes(found_in, metadata, selected_genes, datasets_selected, plot_type, xaxis_type):
    if not found_in.empty:
        ret_str = 'Genes {} were found'.format(selected_genes)

        gene_A, gene_B = found_in["genes"].values[0], found_in["genes"].values[1]

        # create the figure layout
        layout = go.Layout(
            title=" {} vs {} ".format(gene_A, gene_B),
            height=700
        )
        metadata_selected = filter_data(metadata, "Dataset", datasets_selected)

        prcsd_df = found_in.iloc[0]
        prcsd_df = prcsd_df.drop(["genes"]).T.reset_index()
        prcsd_df.columns = ["Sample", "TPM"]
        prcsd_df = sync_df(metadata_selected, prcsd_df)
        prcsd_df.rename(columns={"TPM": gene_A}, inplace=True)
        prcsd_df[gene_B] = found_in[prcsd_df["Sample"].values].iloc[1].values

        hover_data = ["Sample", "Tissue", "Dataset", "subset_name",
                      "NHU_differentiation", "Gender", "TER", "Substrate"]

        config = {"x": gene_A, "y": gene_B, "color": "Dataset",
                  "hover_data": hover_data, "xaxis_type": xaxis_type}
        fig = select_plot_type(prcsd_df, config, plot_type, isComparison=True)

        fig.update_layout(layout)

        if xaxis_type != "log10":
            offset_x = prcsd_df[gene_A].values.max() * 0.015
            offset_y = prcsd_df[gene_B].values.max() * 0.015
            fig.update_xaxes(
                range=[-offset_x, prcsd_df[gene_A].values.max() + offset_x])
            fig.update_yaxes(
                range=[-offset_y, prcsd_df[gene_B].values.max() + offset_y])
        else:
            offset_x = np.log10(list(prcsd_df[gene_A].values)).max() * 0.015
            offset_y = np.log10(list(prcsd_df[gene_B].values)).max() * 0.015
            fig.update_xaxes(
                range=[-offset_x, np.log10(list(prcsd_df[gene_A].values)).max() + offset_x])
            fig.update_yaxes(
                range=[-offset_y, np.log10(list(prcsd_df[gene_B].values)).max() + offset_y])

        return ret_str, fig, prcsd_df

    return "Genes {} were not found in any of the datasets".format(selected_genes), {"data": {}}, found_in


def select_plot_type(df, config, plot_type, isComparison=False, ter=None):
    x, y, color, hover_data, xaxis_type = config["x"], config[
        "y"], config["color"], config["hover_data"], config["xaxis_type"]
    ret_fig = {}
    categorical_order = []

    log = False
    if xaxis_type != "linear":
        log = True
        if not isComparison:
            df.loc[df[df["TPM"] < 1].index.values]["TPM"] = 1.0
        else:
            df.loc[df[df[y] < 1].index.values][y] = 1.0
            df.loc[df[df[x] < 1].index.values][x] = 1.0

    if not isComparison:
        # slightly modified the subtypes and duplicate the ones that have tight TER barrier
        if ter != None and "tight" in ter:
            tight_limit = 500
            ter_col = []
            for _, row in df.iterrows():
                new_name = ""
                if row["isTER"]:
                    if float(row["TER"]) >= tight_limit:
                        new_name = row[x]+"_tight"
                    else:
                        new_name = row[x]+"_non-tight"
                else:
                    new_name = row[x]

                ter_col.append(new_name)

                if new_name not in categorical_order:
                    categorical_order.extend([new_name])

            categorical_order.sort()
            df[x] = ter_col
            color = x

        # do the plot
        if plot_type == "violin":
            ret_fig = px.violin(df, x=x, y=y, color=color,
                                box=True, hover_data=hover_data, log_y=log)
        elif plot_type == "violin_points":
            ret_fig = px.violin(df, x=x, y=y, color=color, box=True,
                                points="all", hover_data=hover_data, log_y=log)
        elif plot_type == "box":
            ret_fig = px.box(df, x=x, y=y, color=color, hover_data=hover_data,
                             log_y=log, category_orders={x: categorical_order})
            ret_fig.update_traces(boxmean="sd")
        elif plot_type == "box_points":
            ret_fig = px.box(df, x=x, y=y, color=color, points="all",
                             hover_data=hover_data, log_y=log,  category_orders={x: categorical_order})
            ret_fig.update_traces(boxmean="sd")
        else:
            ret_fig = px.strip(df, x=x, y=y, color=color, hover_data=hover_data,
                               log_y=log, category_orders={x: categorical_order})
    else:
        ret_fig = px.strip(df, x=x, y=y, color=color,
                           hover_data=hover_data, log_y=log, log_x=log)

    return ret_fig


def change_symbol_ter(df, config, plot_type):
    ter_df = df[df["TER"].astype(np.float16) >= 500]
    if (not ter_df.empty) and (plot_type not in ['violin', 'box']):
        # mark the datasets
        ter_figs = []
        for dataset in ter_df["Dataset"].unique().tolist():
            ter_fig = select_plot_type(
                ter_df[ter_df["Dataset"] == dataset], config, plot_type)
            ter_fig.update_traces(marker_symbol='x-open',
                                  marker_size=13, marker_color="black")
            ter_fig.update_traces(name=ter_fig.data[0].name+"_tight")
            ter_figs.append(ter_fig)
        return ter_figs

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
                ret_string, last_fig, prcsd_data = create_datasets_plot(
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
                ret_string, last_fig, prcsd_data = compare_genes( found_in, metadata, user_input, datasets_selected, plot_type, xaxis_type)

                # normal
                gene_A, gene_B = found_in["genes"].values[0], found_in["genes"].values[1]
                dummy_df = found_in.drop("genes", axis=1).transpose()
                dummy_df.columns = found_in["genes"].values
                pearson_corr, spearman_corr = compute_correlations(dummy_df, gene_A, gene_B)

                # for log10
                dummy_df = pd.DataFrame()
                dummy_df[found_in["genes"].values[0]] = np.log10(list(found_in.iloc[0, 1:].values + 1))
                dummy_df[found_in["genes"].values[1]] = np.log10(list(found_in.iloc[1, 1:].values + 1))
                # we don't need to calculate for spearman, but it's easier to the function
                pearson_corr_log, _ = compute_correlations( dummy_df, gene_A, gene_B, isLog=True)

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
# data_dict = import_data("data/")
# init_callbacks(dash_app, data_dict)

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
            DataTable(id="pearson-table", columns=[], data=[],
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
