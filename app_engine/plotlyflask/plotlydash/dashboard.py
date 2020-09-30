import numpy as np
import pandas as pd
import dash
import dash_table
import dash_html_components as html
import dash_core_components as dcc
from dash.dependencies import Input, Output, State

import plotly.express as px
import plotly.graph_objects as go

from os import path

import time 
from datetime import datetime
import re

def sync_df(df_metadata, df):
    # syncrhonise dfs
    #the df with common samples from both df metadata and the df with tpms
    df_common = filter_data(df, "Sample", list(df_metadata["Sample"].values))
    #for cases when some samples are in a database but not others
    metadata_common = filter_data(df_metadata, "Sample", list(df["Sample"].values)) 

    # sort both by samples so that there are syncd
    metadata_sorted = metadata_common.sort_values("Sample")
    df_common = df_common.sort_values("Sample")
    metadata_sorted["TPM"] = df_common["TPM"].values
    return metadata_sorted

def process_metadata(df_metadata):

    df_metadata = df_metadata.replace("?", np.nan)
    cols = ["NHU_differentiation", "Tissue", "Gender", "Diagnosis", "Substrate"] #cols to replace nan with NA string
    for col in cols:
        df_metadata[col] = df_metadata[col].fillna("Not-known")

    return df_metadata

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
        ret_dict["metadata"] = process_metadata(raw_df)
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
    all_tsv, metadata = pd.DataFrame(), pd.DataFrame()
    ret_dict = {}
    for index, row in master_cv.iterrows():
        filename, description = row["Filename"], row["Description"]
        df = pd.read_csv(base_path+filename, delimiter="\t")
        print("Dataset: {} Shape {}".format(filename, df.shape))
    
        if "metadata" not in filename:
            if not all_tsv.empty:
                all_tsv = pd.merge(all_tsv, df, how="outer")
            else:
                all_tsv = df
        else:
            metadata = df 

    print("Finished loading the data in {}".format(time.time()-start_time))
    ret_dict["metadata"] = process_metadata(metadata)
    ret_dict["all_tsv"] = all_tsv

    start_time = time.time()
    all_tsv.to_csv(base_path + "all_data.tsv",  sep='\t', index = False)
    print("Update .tsv file with the new data in {}".format(time.time() - start_time))
    # return ret_dict
    
def create_datasets_plot(found_in, metadata, user_input, datasets_selected):
    if not found_in.empty:
        ret_str = 'Gene {} was found'.format(user_input)

        # create the figure layout
        layout = go.Layout (
            title = "Swarm plot for {}".format(found_in["genes"].values[0]),
            height= 600
        )
        # prepare datafrarme
        found_in = found_in.drop(["genes"], axis=1).T.reset_index()
        found_in.columns = ["Sample", "TPM"]
        # process metadata and add to the df
        metadata_selected = filter_data(metadata, "Dataset", datasets_selected)
        processed_df = sync_df(metadata_selected, found_in)

        if not datasets_selected:
            return "No datapoints to plot", { "data" : {}}  
        elif len(datasets_selected) == 1:
            # create figure
            tissue_types = []
            for sample in processed_df["Sample"]:
                # tissue_types.append(sample.split("-")[1].strip())
                tissue_type = re.split("\W+", sample)
                print(tissue_type)
                if len(tissue_type) > 1:
                    tissue_types.append(tissue_type[-1])
                else:
                    tissue_types.append(tissue_type[0])

            processed_df["tissue_type"] = tissue_types
            # NOTE: ANDREW  modify this for new label shared_num_same_col
            fig = px.strip(processed_df, x='shared_num_same_col', y='TPM', color="Dataset", 
                    hover_data = ["Tissue", "NHU_differentiation", "Gender", "TER", "Substrate", "Sample", "shared_num_same_col"])
            fig.update_layout(layout)
            return ret_str, fig
        else:
            # create figure
            fig = px.strip(processed_df, x='Dataset', y='TPM', color="Dataset", 
                    hover_data = ["Tissue", "NHU_differentiation", "Gender", "TER", "Substrate"])
            fig.update_layout(layout)
            return ret_str, fig
                
    else:
        return "Gene {} not found in any of the datasets".format(user_input), { "data" : {}} 

def create_medata_plot(found_in, metadata, user_input, datasets_selected,
                            tissue, diagnosis, substrate, ter, gender, xaxis, incl_values):
    if not found_in.empty:
        ret_str = 'Gene {} was found'.format(user_input)
        # create the figure layout
        layout = go.Layout (
            title = "Swarm plot for {}".format(found_in["genes"].values[0]),
            height= 600
        )
        # prepare datafrarme
        found_in = found_in.drop(["genes"], axis=1).T.reset_index()
        found_in.columns = ["Sample", "TPM"]
        # process metadata and add to the df
        metadata_selected = filter_data(metadata, "Dataset", datasets_selected)
        metadata_selected = filter_data(metadata_selected, "Tissue", tissue)
        metadata_selected = filter_data(metadata_selected, "Substrate", substrate)
        metadata_selected = filter_data(metadata_selected, "Gender", gender)

        if ter != None and "tight" in ter:
            metadata_selected = metadata_selected[metadata_selected["TER"].astype(np.float16) >= 500]

        # sort by the incl_cols. The below code will leave only the data that has Y values to the sel columns from incl values
        sel_cols = [col for col in metadata_selected.columns if col in incl_values]
        for sel_col in sel_cols:
            if not metadata_selected.empty:
                metadata_selected = filter_data(metadata_selected, sel_col, ["Y"])

        if metadata_selected.empty:
            return "No TPMs values for gene: {} found with this configurations".format(user_input),  {"data" : {}} 

        # synchronize the metadata df with the tpm dataset if not empty
        metadata_selected = sync_df(metadata_selected, found_in)
    
        # create figure
        fig = px.strip(metadata_selected, x=xaxis, y='TPM', color=xaxis, 
                hover_data = ["Tissue", "NHU_differentiation", "Gender", "TER", "Substrate"])
        fig.update_layout(layout)
        return ret_str, fig
    else:
        return "Gene {} not found in any of the datasets".format(user_input), { "data" : {}}  

def filter_data(df, col, values_to_keep):
    if values_to_keep:
        df = df[df[col].isin(values_to_keep)]
    return df

def create_datasets_tick_box(data_dict):
    # itereate through each of the dictionary and create the html code
    # add calbacks 
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
        html.Br(),
        dcc.RadioItems(
            id = "select-all-none",
            options= [
                {'label': 'All datasets', 'value': 'all'},
                {'label': 'None datasets', 'value': 'none'},
            ]
        ),
    ])

def create_incl_tissue(metadata):
    incl_cols = [col for col in metadata.columns if "incl_" in col]
    options = []
    for incl_col in incl_cols:
        options.append({
            "label": incl_col,
            "value": incl_col
        })
    return html.Div([
        html.H5("Include the following"),
        dcc.Checklist(
            id = "incl-checklist",
            options = options, value = [],
            style = { "column-count": "2"}),
        html.Br()
        ])

def create_gene_search(data_dict):
    # column names should be the same as the ones from the dataframe containing the gen
    return html.Div(id="gene-controls", children = [
                create_datasets_tick_box(data_dict),
                html.Br(),
                html.Button(id='plot-gene', n_clicks=0, children='Plot'),
                html.Div(id='search-result'), 
                html.Br(),
            ])

def metadata_menu(metadata_df):
    # we need to remove some non-numerical values
    cleaned_ter = metadata_df["TER"].dropna().astype(float)
    slider_markers = {}
    for marker in range(0, int(cleaned_ter.max()), 1000):
        slider_markers[marker] = {"label": str(marker)}

    return html.Div(id="metadata-menu", style={"column-count":"1", "margin-top": "40%"},
     children=[
        html.H5("Metadata panel"),
        create_incl_tissue(metadata_df),
        html.Label('Indicate tight barrier? (>500)'),
        dcc.Checklist(id="ter-checklist",
                options=[
                    {'label': 'Tight TER barrier', 'value': 'tight'}
        ]), html.Br(),
        html.Label('Select the X Axis'),
        dcc.Dropdown(
            id="xaxis-dropdown", value="Tissue",
            options=[{"label": i, "value": i} for i in ["NHU_differentiation", "Tissue", "Gender", "Diagnosis", "Substrate", "Dataset"]]
        ), html.Br(),
        html.Label('Tissue'),
        dcc.Dropdown(
            id="tissue-dropdown", value=[metadata_df["Tissue"].dropna().unique()[0]], multi=True,
            options=[{"label": i, "value": i} for i in metadata_df["Tissue"].dropna().unique()]
        ),html.Br(),
        html.Label('Diagnosis'),
        dcc.Dropdown(
            id="diagnosis-dropdown", value=[metadata_df["Diagnosis"].dropna().unique()[0]], multi=True,
            options=[{"label": i, "value": i} for i in metadata_df["Diagnosis"].dropna().unique()]
        ), html.Br(),
        html.Label('Substrate'),
        dcc.Dropdown(
            id="substrate-dropdown", value=[metadata_df["Substrate"].dropna().unique()[0]], multi=True,
            options=[{"label": i, "value": i} for i in metadata_df["Substrate"].dropna().unique()]
        ), html.Br(),
        html.Label('Gender'),
        dcc.Dropdown(
            id="gender-radiobtn", value=metadata_df["Gender"].dropna().unique(), multi=True,
            options= [{"label": i, "value": i} for i in metadata_df["Gender"].dropna().unique()]
        ), html.Br(),
        html.Label("After you selected the configuration from above press the below button to plot"), 
        html.Br(),
        html.Button(id='metadata-plot', n_clicks=0, children='Plot'),
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
         [Output('search-result', "children"), Output('dataset-plots', "figure")],
         [Input('plot-gene', 'n_clicks'), Input('metadata-plot', 'n_clicks'),
          Input("gene-input","n_submit")],
         [State("gene-input", "value"),  State("data-checklist", "value"),
         State("tissue-dropdown","value"), State("diagnosis-dropdown","value"),
         State("substrate-dropdown","value"), State("ter-checklist","value"), State("gender-radiobtn","value"), State("xaxis-dropdown", "value"), State("incl-checklist", "value")]
    )
    def button_router(btn_1, btn_2, n_submit, user_input, datasets_selected, 
                        tissue, diagnosis, substrate, ter, gender, xaxis, incl_values):
        ctx = dash.callback_context

        # get the button_id
        if not ctx.triggered:
            button_id = 'No clicks yet'
        else:
            button_id = ctx.triggered[0]['prop_id'].split('.')[0]

        prcsd_str = user_input.strip()
        if prcsd_str:
            all_tsv = data_dict["all_tsv"]
            metadata = data_dict["metadata"]
            found_in = all_tsv[all_tsv["genes"].str.fullmatch("\\b"+prcsd_str+"\\b", case = False, na=False)]

            if button_id == "plot-gene" or button_id == "gene-input":
                print("Dataset plot has been trieggered")  
                ret_string, last_fig = create_datasets_plot(found_in, metadata, prcsd_str, datasets_selected)
                return ret_string, last_fig
            else:
                print("Metadata plot has been triggered")
                ret_string, last_fig = create_medata_plot(found_in, metadata, prcsd_str, datasets_selected, tissue, diagnosis, substrate, ter, gender, xaxis, incl_values)
                return ret_string, last_fig
        return "Search for gene and click submit", { "data" : {} }

    @dash_app.callback(
        Output("data-checklist", "value"),
        [Input("select-all-none", "value")]
    )
    def select_all_none_dasets(radio_value):    
        if radio_value == "none":
            return []
        else:
            return list(data_dict["metadata"]["Dataset"].unique())

    @dash_app.callback(
        Output("output-export-csv", 'children'),
        [Input("export-plot", "n_clicks")],
        [State("plot-width", "value"), State("plot-height", "value"), State("plot-scale", "value"), State("dataset-plots", "figure")]
    )
    def export_plot(btn_1, width, height, scale, figure):
        path =  "./exported_data/"
        ret_string = "There is no figure to save"
        if figure["data"]:
            fig = go.Figure(figure)
            name = "{}jbu-viz-{}{}".format(path, datetime.now().strftime('%d-%m-%Y--%H:%M:%S'), ".pdf")
            fig.write_image(name, width=width, height=height, scale=scale) 
            ret_string = "Plot saved to {}".format(name)
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
            export_panel()
        ]),
    ])
    return dash_app.server




