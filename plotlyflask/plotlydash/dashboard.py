import numpy as np
import pandas as pd
import dash
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
        ret_dict["all_tsv"] = pd.read_csv(base_path + "all_data.tsv", delimiter="\t")
        raw_df = pd.read_csv(base_path + "00_JBU_sequencing_metadata.tsv", delimiter="\t")
        # we need to do some pre-processing to duplicate the data _incl for fields
        ret_dict["metadata"], ret_dict["all_tsv"] = process_data(raw_df, ret_dict["all_tsv"])
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
    all_tsv.to_csv(base_path + "all_data.tsv",  sep='\t', index = False)
    print("Update .tsv file with the new data in {}".format(time.time() - start_time))

def process_data(df_metadata, all_tsv):
    """ Process the metadata and the all_tsv file. The are two main operations happening:
    1. Replaing in specified cols the nan values and "?" with Not-known
    2. Fakes the _incl columns with values of "Y" to be as datasets. This involved to duplicate both the samples in metadata and in the `all_data.tsv`. We've marked the samples that are duplicated by the following rule "_incl_" + "including column initials". 
    Args:
        df_metadata ([DataFrame]): The medata
        all_tsv ([DataFrame]): TPM dataframe
    Returns:
        (DataFrame, DataFrame): The two dataframes updated
    """

    df_metadata = df_metadata.replace("?", np.nan)
    #cols to replace nan with NA string
    cols = ["NHU_differentiation", "Tissue", "Gender", "Diagnosis", "Substrate"]
    for col in cols:
        df_metadata[col] = df_metadata[col].fillna("Not-known")

    #add the incl columns to the datasets and duplicate the samples
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
        duplicated_samples = [sample+initials for sample in metadata_selected["Sample"].values]
        metadata_selected = metadata_selected.drop("Sample", axis=1)
        metadata_selected["Sample"] = duplicated_samples
        metadata_selected["Dataset"] = incl_col
        metadata_selected["subset_name"] = incl_col
        new_rows_metadata = new_rows_metadata.append(metadata_selected, ignore_index=True)

        # update all_tsv
        all_tsv_selected = all_tsv[selected_samples]
        all_tsv_selected.columns = duplicated_samples
        new_cols_all_tsv = pd.concat([new_cols_all_tsv, all_tsv_selected], axis=1)
        
    df_metadata = pd.merge(df_metadata, new_rows_metadata, how="outer")
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
    #the df with common samples from both df metadata and the df with tpms
    df_common = filter_data(df, "Sample", list(df_metadata["Sample"].values))
    
    #for cases when some samples are in a database but not others
    metadata_common = filter_data(df_metadata, "Sample", list(df["Sample"].values)) 

    # sort both by samples so that there are syncd
    metadata_sorted = metadata_common.sort_values("Sample")
    df_common = df_common.sort_values("Sample")
    metadata_sorted["TPM"] = df_common["TPM"].values
    return metadata_sorted

def samples_to_remove(df):
    #used to remove _incl samples when exporting the .csv files
    ret_samples = []
    # samples duplicated by tpm
    duplicated_samples = df[df["TPM"].duplicated()]["Sample"].values
    # but we need to check if these are the duplicates of the _incl (maybe in the future there will be samples that have the exact TPM value)
    ret_samples = [sample for sample in duplicated_samples if "_incl_" in sample]
    return ret_samples

##### Plotting functions #####
def create_datasets_plot(found_in, metadata, user_input, datasets_selected, ter, plot_type, xaxis, xaxis_type):
    if not found_in.empty:
        ret_str = 'Gene {} was found'.format(user_input)

        # create the figure layout
        layout = go.Layout (
            title = "Swarm plot for {}".format(found_in["genes"].values[0]),
            height= 700
        )
        # prepare datafrarme
        found_in = found_in.drop(["genes"], axis=1).T.reset_index()
        found_in.columns = ["Sample", "TPM"]

        # process metadata and add to the df
        metadata_selected = filter_data(metadata, "Dataset", datasets_selected)
        processed_df = sync_df(metadata_selected, found_in)

        hover_data = ["Sample",  "shared_num_same_col", "Tissue", "Dataset","NHU_differentiation", "Gender", "TER", "Substrate"]
        if not processed_df.empty:
            if not datasets_selected:
                return "Select a dataset first", { "data" : {}}, found_in

            # create the main trace
            
            config={"x":xaxis , "y":"TPM", "color":"Dataset", "hover_data": hover_data, "xaxis_type": xaxis_type}
            fig = select_plot_type(processed_df, config, plot_type)

            # add the ter symbols if needed 
            if ter != None and "tight" in ter:
                ter_figs = change_symbol_ter(processed_df, config, plot_type)
                if ter_figs is not None:
                    for ter_fig in ter_figs:
                        fig.add_trace(ter_fig.data[0])

            # if the maximum value for tpm is 25 change the tick
            if xaxis_type == "linear":
                if processed_df["TPM"].values.max() <= 25:
                    fig.update_yaxes(range=[0, 25])
                else:
                    fig.update_yaxes(range=[0, processed_df["TPM"].values.max()])

            fig.update_layout(layout)
            return ret_str, fig, processed_df

    return "Gene {} not found in any of the datasets".format(user_input), { "data" : {}}, found_in 

def select_plot_type(df, config, plot_type):
    x, y, color, hover_data, xaxis_type = config["x"], config["y"], config["color"], config["hover_data"], config["xaxis_type"]
    ret_fig = {}

    log = False
    if xaxis_type != "linear":
        log = True
        df.loc[df[df["TPM"] < 1].index.values]["TPM"] = 1.0

    if plot_type == "violin":
        ret_fig = px.violin(df, x=x, y=y, color=color, box=True, hover_data = hover_data, log_y=log)
    elif plot_type == "violin_points":
        ret_fig = px.violin(df, x=x, y=y, color=color, box=True,  points="all", hover_data = hover_data, log_y=log)
    elif plot_type == "box":
        ret_fig = px.box(df, x=x, y=y, color=color, hover_data = hover_data, log_y=log)
        ret_fig.update_traces(boxmean="sd")
    elif plot_type == "box_points":
        ret_fig = px.box(df, x=x, y=y, color=color, points="all", hover_data = hover_data, log_y=log)
        ret_fig.update_traces(boxmean="sd")
    else:
        ret_fig = px.strip(df, x=x, y=y, color=color, hover_data = hover_data, log_y=log)
        
    return ret_fig

def change_symbol_ter(df, config, plot_type):
    ter_df = df[df["TER"].astype(np.float16) >= 500]
    if (not ter_df.empty) and (plot_type not in ['violin', 'box']):
        # mark the datasets
        ter_figs = []
        for dataset in ter_df["Dataset"].unique().tolist():
            ter_fig = select_plot_type(ter_df[ter_df["Dataset"] == dataset], config, plot_type)
            ter_fig.update_traces(marker_symbol='x-open', marker_size=13, marker_color="black")
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
        html.H5("Change the value in the text box to see callbacks in action!"),
        dcc.Checklist(
            id = "data-checklist",
            options = options, value = datasets,
            style = { "column-count": "2"}),
    ])

#TODO: Rename this to create dataset panel
def create_gene_search(data_dict):
    # column names should be the same as the ones from the dataframe containing the gen
    return html.Div(id="gene-controls", children = [
                create_datasets_tick_box(data_dict),
                html.Br(),
                dcc.RadioItems(
                    id = "select-all-none",
                    options= [
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

def metadata_menu(metadata_df):
    # we need to remove some non-numerical values
    x_axis_options = ["NHU_differentiation", "Tissue", "Gender", "Diagnosis", "Substrate", "Dataset", "subset_name"]
    return html.Div(id="metadata-menu", style={"column-count":"1", "margin-top": "40%"},
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
            options = [
                {'label':'Swarm', 'value':'swarm',},
                {'label':'Violin', 'value':'violin'},
                {'label': 'Vilion and points', 'value': 'violin_points'},
                {'label':'Box', 'value':'box'},
                {'label':'Box and points', 'value':'box_points'}
            ]
        ),
        html.H6('Select the X Axis type'),
        dcc.RadioItems(
            id="select-xaxis-type", value='linear',
            options = [
                {'label':'Linear', 'value':'linear'},
                {'label':'Log10 ', 'value':'log10'}
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
        html.Div( style = { "column-count": "2"}, children=[
            html.Button(id='export-plot', n_clicks=0, children='Export Plot'),
            html.Button(id='export-data', n_clicks=0, children='Export CSV')
            ]),
        html.Div(id="output-export-csv")
    ])

def init_callbacks(dash_app, data_dict):
    # latest fig to save
    @dash_app.callback(
         [Output('search-result', "children"), Output('dataset-plots', "figure"), Output('intermediate-value', 'children')],
         [Input('plot-gene', 'n_clicks'),
          Input("gene-input","n_submit")],
         [State("gene-input", "value"),  State("data-checklist", "value"), State("ter-checklist","value"), State("xaxis-dropdown", "value"), State("select-plot-type","value"), State("select-xaxis-type","value") ]
    )
    def button_router(btn_1, n_submit, user_input, datasets_selected, ter, xaxis, plot_type, xaxis_type):

        ret_data = pd.DataFrame().to_json(date_format='iso', orient='split') 
        prcsd_str = user_input.strip()
        if prcsd_str:
            all_tsv, metadata = data_dict["all_tsv"], data_dict["metadata"]
            # search for the gene
            found_in = all_tsv[all_tsv["genes"].str.fullmatch("\\b"+prcsd_str+"\\b", case = False, na=False)]
            # create the figure
            ret_string, last_fig, prcsd_data = create_datasets_plot(found_in, metadata, prcsd_str, datasets_selected, ter, plot_type, xaxis, xaxis_type)
            # transformed the filtered data to json
            ret_data = prcsd_data.to_json(date_format='iso', orient='split')

            return ret_string, last_fig, ret_data
        return "Search for gene and click submit", { "data" : {} }, ret_data

    @dash_app.callback(
        Output("data-checklist", "value"),
        [Input("select-all-none", "value")],
        [State("data-checklist", "value")]
    )
    def select_all_none_dasets(radio_value, datasets):    
        if radio_value == "none":
            return []
        elif radio_value == "excl":     
            incl_cols = [col for col in data_dict["metadata"].columns if "incl_" in col]   
            for incl_col in incl_cols:
                if incl_col in datasets:
                    datasets.remove(incl_col)
            return datasets
        elif radio_value == "incl":
            incl_cols = [col for col in data_dict["metadata"].columns if "incl_" in col]
            # include only the incl_cols that were not selected
            datasets+=list(set(incl_cols) - set(datasets))
            return datasets
        else:
            return list(data_dict["metadata"]["Dataset"].unique())

    @dash_app.callback(
        Output("output-export-csv", 'children'),
        [Input("export-plot", "n_clicks"), Input("export-data", "n_clicks"), Input('intermediate-value', 'children')],
        [State("plot-width", "value"), State("plot-height", "value"), State("plot-scale", "value"), State("dataset-plots", "figure"), State("data-checklist", "value")]
    )
    def export_plot(btn_1, btn_2, data_json, width, height, scale, figure, datasets):
        # checking os
        if name == 'nt':
            path =  "exported_data\\"
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
                filepath = "{}plot-jbu-viz-{}{}".format(path, datetime.now().strftime('%d-%m-%Y--%H-%M-%S'), ".pdf")
                fig.write_image(filepath, width=width, height=height, scale=scale) 
                ret_string = "Plot saved to {}".format(filepath)
        elif button_id == "export-data":
            figure_data_df = pd.read_json(data_json, orient='split')
            filepath = "{}data-jbu-viz-{}{}".format(path, datetime.now().strftime('%d-%m-%Y--%H-%M-%S'), ".csv")

            # we need to remove the duplicated samples
            unique_subsets = figure_data_df["subset_name"].unique()
            export_df = pd.DataFrame()
            for subset_name in unique_subsets:
                dummy_df = figure_data_df[figure_data_df["subset_name"] == subset_name]["TPM"].T
                # we need to rename the series so that has a unique name, otherwise the below line won't work
                export_df = pd.concat([export_df, dummy_df.rename(subset_name)], axis=1)

            # we need to put the nan values at the bottom and to that we need to sort each column individualy and re-concat
            export_df = pd.concat([export_df[col].sort_values().reset_index(drop=True) for col in export_df], axis=1, ignore_index=True)
            export_df.columns = [subset+"_TPM" for subset in unique_subsets] 

            export_df.to_csv(filepath, index=False)
            ret_string = "Data saved to {}".format(filepath)
        return ret_string
                
def init_dashboard(server):
    """Create a Plotly Dash dashboard."""
    external_stylesheets = ['https://codepen.io/chriddyp/pen/bWLwgP.css']
    dash_app = dash.Dash(
        server=server,
        routes_pathname_prefix='/dashapp/',
        external_stylesheets=external_stylesheets
    )

    # update_all_datasets("data/")
    data_dict = import_data("data/")
    init_callbacks(dash_app, data_dict)

    # Create Layout
    dash_app.layout = html.Div(children=[
        html.H1(children='JBU visualisation tool'),
        html.Div(id="parrent-all", children=[
            html.Div(id="gene-search", children=[
                html.H5("Enter the gene you want to analyse"),
                html.Div(["Gene name: ",
                    dcc.Input(id='gene-input', value="", type='text')]),
                html.Br(),
            ]),
            html.Div(id="parent-controllers", style = { "column-count": "2", "width":"fit-content"}, 
            children=[
                create_gene_search(data_dict), 
                metadata_menu(data_dict["metadata"]),
            ]),
            dcc.Graph(
                id='dataset-plots',
                figure = { "data" : {},
                "layout": { "title": "No Gene Search", "height": 300 }
                }),
            export_panel(),
            html.Div(id='intermediate-value', style={'display': 'none'}) #used to store data between callbacks
        ]),
    ])
    return dash_app.server

