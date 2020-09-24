import numpy as np
import pandas as pd
import dash
import dash_table
import dash_html_components as html
import dash_core_components as dcc
from dash.dependencies import Input, Output, State

import plotly.express as px
import plotly.graph_objects as go

import time 

def import_files(base_path):

    # open the master .csv file, containing the contents of the folder
    master_cv = pd.read_csv(base_path + "master.tsv", delimiter="\t")

    # iterate through the csv file and create a list of dictionaries with the following: <key> - name of the file, <value> object which contains a dataframe and description from the .csv file.
    ret_dict = {}
    all_tsv = pd.DataFrame()
    combined_genes = pd.DataFrame() 
    for index, row in master_cv.iterrows():
        filename, description = row["Filename"], row["Description"]
        df = pd.read_csv(base_path+filename, delimiter="\t")
        print("Dataset: {} Shape {}".format(filename, df.shape))
        ret_dict[filename] = {
            "df": df,
            "description": description 
            }
    
        # df["dataset"] = filename.replace(".tsv", "")
        if "metadata" not in filename:
            if not all_tsv.empty:
                all_tsv = pd.merge(all_tsv, df, how="outer")
            else:
                all_tsv = df

        # get gene column
        if "genes" in df.columns.values:
            genes = pd.DataFrame(df["genes"])
            genes["dataset"] = filename
            combined_genes = combined_genes.append(genes)

    ret_dict["all_tsv"] = {
        "df": all_tsv,
        "description": "all tsv files in one place"
    }
    return ret_dict, combined_genes
    
def create_swarm_figs(metadata, found_in, datasets_selected):
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
    return fig

def create_dfs_tick_box(data_dict):
    # itereate through each of the dictionary and create the html code
    # add calbacks 
    options = []
    datasets = data_dict["00_JBU_sequencing_metadata.tsv"]["df"]["Dataset"].unique()
    for dataset in datasets:
        options.append({
            "label": dataset,
            "value": dataset,
        })
    return html.Div([
        html.H6("Change the value in the text box to see callbacks in action!"),
        dcc.Checklist(
            id = "data-checklist",
            options = options,
            value = datasets,
            style = { "column-count": "2", "width": "50%"}),
        html.Br(),
    ])

def create_gene_search(data_dict):
    # column names should be the same as the ones from the dataframe containing the gen
    return html.Div(id="gene-controls", children = [
                html.H6("Find a gene in the dataset"),
                html.Div([
                    "Gene name: ",
                    dcc.Input(id='gene-input', value="", type='text')]),
                html.Br(),
                create_dfs_tick_box(data_dict),
                html.Br(),
                html.Button(id='plot-gene', n_clicks=0, children='Plot'),
                html.Div(id='search-result')
            ])

def metadata_menu(metadata_df):
    return html.Div(id="metadata-menu", children=[
        html.H6("Select sample(s)"),
        html.Label('Select Sample'),
        dcc.Dropdown(
            id="sample-dropdown",
            options=[{"label": i, "value": i} for i in metadata_df["Sample"]],
            value="Select samples",
            multi=True
        ),
        html.Br(),
        html.Label('Choose Tissue'),
        dcc.Dropdown(
            id="tissue-dropdown",
            options=[{"label": i, "value": i} for i in metadata_df["Tissue"].drop_duplicates()],
            value="Select samples",
            multi=True
        ),
        html.Br(),
        html.Label('Select Diagnosis'),
        dcc.Dropdown(
            id="diagnosis-dropdown",
            options=[{"label": i, "value": i} for i in metadata_df["Diagnosis"].drop_duplicates()],
            value="Select samples",
            multi=True
        ),
        html.Br(),
        html.Label('Choose Substrate'),
        dcc.Dropdown(
            id="substrate-dropdown",
            options=[{"label": i, "value": i} for i in metadata_df["Substrate"].drop_duplicates()],
            value="Select samples",
            multi=True
        ),
        html.Br(),
        html.Label('Choose Ter'),
        dcc.Dropdown(
            id="ter-dropdown",
            options=[{"label": i, "value": i} for i in metadata_df["TER"].drop_duplicates()],
            value="Select samples",
            multi=True
        ),
        html.Br(),
        html.Label('Choose Gender'),
        dcc.Checklist(
            id="gender-checklist",
            options= [
                {'label': 'Male', 'value': 'male'},
                {'label': 'Female', 'value': 'female'},
                {'label': 'Other', 'value': 'other'}
            ]
        )
        # dropdown for substrate
        # dropdown fro ter
        # checkbox for gender
    ])
        
def init_callbacks(dash_app, data_dict, combined_genes):
    # @dash_app.callback(
    # Output('plot-gene', "children"),
    # [])
    # def show_selected_data(input_value):
    #     return "Re-Plot"

    @dash_app.callback(
        [Output('search-result', "children"),
        Output('dataset-plots', "figure"),
        Output('plot-gene', "children")],
        [Input('plot-gene', 'n_clicks')],
        [State("gene-input", "value"),
        State("data-checklist", "value")])
    def show_gene_plots(n_clicks, user_input, datasets_selected):
        if user_input.strip():
            all_tsv = data_dict["all_tsv"]["df"]

            start_time = time.time()
            found_in = all_tsv[all_tsv["genes"].str.match("\\b"+user_input+"\\b", case = False, na=False)]
            print("Finished search in all tsv in: {:.4f}".format(time.time()-start_time))
           
            start_time = time.time()
            found_in_df = combined_genes[combined_genes["genes"].str.match("\\b"+user_input+"\\b", case = False, na=False)]
            print("Finished search in genes seperate array in: {:.4f}".format(time.time()-start_time))
            datasets = " ".join(found_in_df["dataset"].astype(str))
            if not found_in.empty:
                ret_str = 'Gene {} found in dataset(s): {}'.format(user_input, datasets)

                fig = create_swarm_figs(data_dict["00_JBU_sequencing_metadata.tsv"]["df"], found_in, datasets_selected)
                return ret_str, fig, "Plot"
            else:
                return "Gene {} not found in any of the datasets".format(user_input), {}
        return "Search for gene and click submit", {}, "Plot" 

def init_dashboard(server):
    """Create a Plotly Dash dashboard."""
    external_stylesheets = ['https://codepen.io/chriddyp/pen/bWLwgP.css']
    dash_app = dash.Dash(
        server=server,
        routes_pathname_prefix='/dashapp/',
        external_stylesheets=external_stylesheets
    )
    data_dict, combined_genes = import_files("/Users/vlad/Documents/Code/York/visualisation/data/")
    init_callbacks(dash_app, data_dict, combined_genes)

    # Create Layout
    dash_app.layout = html.Div(children=[
        html.H1(children='JBU visualisation tool'),
        html.Div(id="parrent-all", children=[
            html.Div(id="parent-controllers", style = { "column-count": "2", "width":"fit-content"}, 
            children=[
                create_gene_search(data_dict), 
                metadata_menu(data_dict["00_JBU_sequencing_metadata.tsv"]["df"]),
            ]),
            dcc.Graph(
                id='dataset-plots',
                figure = { "data" : {},
                "layout": { "title": "No Gene Search", "height": 300 }
                })
        ]),
    ])
    return dash_app.server

def create_test_plot():

    df = pd.DataFrame({
    "Fruit": ["Apples", "Oranges", "Bananas", "Apples", "Oranges", "Bananas"],
    "Amount": [4, 1, 2, 2, 4, 5],
    "City": ["SF", "SF", "SF", "Montreal", "Montreal", "Montreal"]
    })

    fig = px.bar(df, x="Fruit", y="Amount", color="City", barmode="group")

    return dcc.Graph(
            id='example-graph',
            figure=fig
    )






