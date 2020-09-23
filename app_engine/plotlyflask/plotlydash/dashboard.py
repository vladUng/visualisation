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
    # all_tsv = data_dict["all_tsv"]["df"]
    # metadata = data_dict["00_JBU_sequencing_metadata.tsv"]["df"]
    layout = go.Layout (
        title = "Swarm plot for {}".format(found_in["genes"].values[0]),
        height= 600
        # width = 700,
    )

    found_in = found_in.drop(["genes"], axis=1).T.reset_index()
    found_in.columns = ["Sample", "TPM"]

    metadata_selected = metadata[metadata["Sample"].isin(found_in["Sample"])]
    procesed_df = pd.concat([found_in, metadata_selected], axis=1)

    # drop_datasets = list(set(datasets_selected) & set(procesed_df["Dataset"].values))

    procesed_df = procesed_df[procesed_df["Dataset"].isin(datasets_selected)]
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
    print(datasets)
    return html.Div([
        html.H6("Change the value in the text box to see callbacks in action!"),
        dcc.Checklist(
            id = "data-checklist",
            options = options,
            value = datasets,
            style = {
               "column-count": "2"
            },
        ),
        html.Br(),
        html.Div(id="data-selected")
    ])

def create_gene_search(data_dict):
    # column names should be the same as the ones from the dataframe containing the gen
    return html.Div([
                html.H6("Find a gene in the dataset"),
                html.Div(["Gene name: ",
                        dcc.Input(id='gene-input', value="", type='text'),
                        html.Button(id='search-button-gene', n_clicks=0, children='Submit')],
                        style = {"column-count": "2"}
                        ),
                html.Br(),
                create_dfs_tick_box(data_dict),
                html.Br(),
                html.Div(id='search-result'),
                html.Br(),
                dcc.Graph(
                    id='dataset-plots',
                    figure = { "data" : {},
                    "layout": {
                        "title": "My Dash Graph",
                        "height": 400,  # px
                    }})
            ])
            # dash_table.DataTable(
            #     id='found-dataset',
            #     columns=[{"name": i, "id": i} for i in df.columns],
            #     data=df.to_dict('records'),
            #     style_cell={'textAlign': 'center','min-width':'50px', 'width':'200px'})

def metadata_menu(metadata_df):
    return html.Div(children=[
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
    @dash_app.callback(
        Output(component_id='data-selected', component_property="children"),
        [Input(component_id='data-checklist', component_property='value')])
    def show_selected_data(input_value):
        return 'Output: {}'.format(input_value)

    # @dash_app.callback(
    #     Output(component_id='found-dataset', component_property="data"),
    #     [Input(component_id='gene-input', component_property='value')])
    def search_gene(user_input):
        genes_found, datasets = "", ""
        found_in_df = pd.DataFrame()
        if (user_input.strip()): 
            found_in_df = combined_genes[combined_genes["genes"].str.contains(user_input, case = False, na=False)]
            genes_found = " ".join(found_in_df["genes"].astype(str))
            datasets = " ".join(found_in_df["dataset"].astype(str))
        # return 'Gene(s) found: {} in dataset(s): {}'.format(genes_found, datasets)
        return found_in_df.to_dict('records')
        
    @dash_app.callback(
        [Output('search-result', "children"),
        Output('dataset-plots', "figure")],
        [Input('search-button-gene', 'n_clicks')],
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
                return ret_str, fig
            else:
                return "Gene {} not found in any of the datasets".format(user_input), {}
        return "Search for gene and click submit", {} 

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
        html.H1(children='Hello Dash'),
        html.Div(children='''
            Dash: A web application framework for Python.
        '''),
        html.Div(children=[
        create_gene_search(data_dict), 
        ]),
        metadata_menu(data_dict["00_JBU_sequencing_metadata.tsv"]["df"]),
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






