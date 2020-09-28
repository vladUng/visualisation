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
        metadata_selected = metadata[metadata["Sample"].isin(found_in["Sample"])]
        procesed_df = pd.concat([found_in, metadata_selected], axis=1)
        procesed_df = procesed_df[procesed_df["Dataset"].isin(datasets_selected)]
        # create figure
        fig = px.strip(procesed_df, x='Dataset', y='TPM', color="Dataset", 
                hover_data = ["Tissue", "NHU_differentiation", "Gender", "TER", "Substrate"])
        fig.update_layout(layout)
        return ret_str, fig
    else:
        return "Gene {} not found in any of the datasets".format(user_input), {}   

def create_medata_plot(found_in, metadata, user_input, datasets_selected,
                            tissue, diagnosis, substrate, ter, gender, xaxis):
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

        # synchronize the metadata df with the tpm dataset
        found_in = filter_data(found_in, "Sample", list(metadata_selected["Sample"].values))
        metadata_selected= metadata_selected.sort_values("Sample")
        found_in = found_in.sort_values("Sample")
        metadata_selected["TPM"] = found_in["TPM"]

        # create figure
        fig = px.strip(metadata_selected, x=xaxis, y='TPM', color=xaxis, 
                hover_data = ["Tissue", "NHU_differentiation", "Gender", "TER", "Substrate"])
        fig.update_layout(layout)
        return ret_str, fig
    else:
        return "Gene {} not found in any of the datasets".format(user_input), {}   

def filter_data(df, col, values_to_keep):
    if values_to_keep:
        df = df[df[col].isin(values_to_keep)]
    return df

def create_dfs_tick_box(data_dict):
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
        html.H4("Change the value in the text box to see callbacks in action!"),
        dcc.Checklist(
            id = "data-checklist",
            options = options,
            value = datasets,
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

def create_gene_search(data_dict):
    # column names should be the same as the ones from the dataframe containing the gen
    return html.Div(id="gene-controls", children = [
                html.H3("Find a gene in the dataset"),
                html.Div([
                    "Gene name: ",
                    dcc.Input(id='gene-input', value="", type='text')]),
                html.Br(),
                create_dfs_tick_box(data_dict),
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
        html.Label('Choose Tissue'),
        dcc.Dropdown(
            id="tissue-dropdown", value=[metadata_df["Tissue"].dropna().unique()[0]], multi=True,
            options=[{"label": i, "value": i} for i in metadata_df["Tissue"].dropna().unique()]
        ),
        html.Br(),
        html.Label('Select Diagnosis'),
        dcc.Dropdown(
            id="diagnosis-dropdown", value=[metadata_df["Diagnosis"].dropna().unique()[0]], multi=True,
            options=[{"label": i, "value": i} for i in metadata_df["Diagnosis"].dropna().unique()]
        ),
        html.Br(),
        html.Label('Choose Substrate'),
        dcc.Dropdown(
            id="substrate-dropdown", value=[metadata_df["Substrate"].dropna().unique()[0]], multi=True,
            options=[{"label": i, "value": i} for i in metadata_df["Substrate"].dropna().unique()]
        ),
        html.Br(),
        html.Label('Choose Ter'),
        html.Div([
            dcc.RangeSlider(
                id="ter-range",
                min = cleaned_ter.min(),
                max = cleaned_ter.max(),
                step = 100,
                marks = slider_markers,
                value=[cleaned_ter.median() - 400, cleaned_ter.median() + 400], 
                allowCross = False
            ),
            html.Div(id='output-container-range-slider')
        ]),
        html.Br(),
        html.Label('Choose Gender'),
        dcc.Dropdown(
            id="gender-radiobtn", value=metadata_df["Gender"].dropna().unique(), multi=True,
            options= [{"label": i, "value": i} for i in metadata_df["Gender"].dropna().unique()]
        ),
        html.Label('Choose the X Axis'),
        dcc.Dropdown(
            id="xaxis-dropdown", value="Tissue",
            options=[{"label": i, "value": i} for i in ["NHU_differentiation", "Tissue", "Gender", "Diagnosis", "Substrate", "Dataset"]]
        ), 
        html.Label('Curate the list'),
        dcc.Checklist(
            id="curate-list",
            options= [
                {'label': 'All datasets', 'value': 'all'},
                {'label': 'None datasets', 'value': 'none'},
            ]
        ),
        html.Label('Test'),
        dcc.RadioItems(
            id = "test-radio",
            options= [
                {'label': 'All datasets', 'value': 'all'},
                {'label': 'None datasets', 'value': 'none'},
            ]
        ),
        # html.Br(),
        # html.Label('Choose plot colouring'),
        # dcc.Dropdown(
        #     id="colouring-dropdown", value="Dataset",
        #     options=[{"label": i, "value": i} for i in ["NHU_differentiation", "Tissue", "Gender", "Diagnosis", "Substrate", "Dataset"]]
        # ), html.Br(),
        html.Label("After you selected the configuration from above press the below button to plot"), 
        html.Br(),
        html.Button(id='metadata-plot', n_clicks=0, children='Plot'),
    ])
        
def init_callbacks(dash_app, data_dict):
    @dash_app.callback(
         [Output('search-result', "children"), Output('dataset-plots', "figure")],
         [Input('plot-gene', 'n_clicks'), Input('metadata-plot', 'n_clicks')],
         [State("gene-input", "value"),  State("data-checklist", "value"),
         State("tissue-dropdown","value"), State("diagnosis-dropdown","value"),
         State("substrate-dropdown","value"), State("ter-range","value"), State("gender-radiobtn","value"), State("xaxis-dropdown", "value")]
    )
    def button_router(btn_1, btn_2, user_input, datasets_selected, 
                        tissue, diagnosis, substrate, ter, gender, xaxis):
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
            found_in = all_tsv[all_tsv["genes"].str.match("\\b"+prcsd_str+"\\b", case = False, na=False)]

            if button_id == "plot-gene":
                print("Dataset plot has been trieggered")  
                ret_string, fig = create_datasets_plot(found_in, metadata, prcsd_str, datasets_selected)
                return ret_string, fig
            else:
                print("Metadata plot has been triggered")
                ret_string, fig = create_medata_plot(found_in, metadata, prcsd_str, datasets_selected, 
                                                        tissue, diagnosis, substrate, ter, gender, xaxis)
                return ret_string, fig

        return "Search for gene and click submit", {}

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
    Output('output-container-range-slider', 'children'),
    [Input('ter-range', 'value')])
    def update_output(value):
        return 'You have selected the following Range "{}"'.format(value)

def init_dashboard(server):
    """Create a Plotly Dash dashboard."""
    external_stylesheets = ['https://codepen.io/chriddyp/pen/bWLwgP.css']
    dash_app = dash.Dash(
        server=server,
        routes_pathname_prefix='/dashapp/',
        external_stylesheets=external_stylesheets
    )

    # update_all_datasets("/Users/vlad/Documents/Code/York/visualisation/data/")
    data_dict = import_data("/Users/vlad/Documents/Code/York/visualisation/data/")
    init_callbacks(dash_app, data_dict)

    # Create Layout
    dash_app.layout = html.Div(children=[
        html.H1(children='JBU visualisation tool'),
        html.Div(id="parrent-all", children=[
            html.Div(id="parent-controllers", style = { "column-count": "2", "width":"fit-content"}, 
            children=[
                create_gene_search(data_dict), 
                metadata_menu(data_dict["metadata"]),
            ]),
            dcc.Graph(
                id='dataset-plots',
                figure = { "data" : {},
                "layout": { "title": "No Gene Search", "height": 300 }
                })
        ]),
    ])
    return dash_app.server




