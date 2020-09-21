import numpy as np
import pandas as pd
import dash
import dash_table
import dash_html_components as html
import dash_core_components as dcc
from dash.dependencies import Input, Output

import plotly.express as px

def import_files(base_path):

    # open the master .csv file, containing the contents of the folder
    master_cv = pd.read_csv(base_path + "master.tsv", delimiter="\t")

    # iterate through the csv file and create a list of dictionaries with the following: <key> - name of the file, <value> object which contains a dataframe and description from the .csv file.
    ret_dict = {}
    combined_genes = pd.DataFrame() 
    for _, row in master_cv.iterrows():
        filename, description = row["Filename"], row["Description"]
        df = pd.read_csv(base_path+filename, delimiter="\t")
        ret_dict[filename] = {
            "df": df,
            "description": description 
            }

        # get gene column
        if "genes" in df.columns.values:
            genes = pd.DataFrame(df["genes"])
            genes["dataset"] = filename
            combined_genes = combined_genes.append(genes)
    
    return ret_dict, combined_genes
    
def create_dfs_tick_box(data_dict):
    # itereate through each of the dictionary and create the html code
    # add calbacks 
    options = []
    for key, _ in data_dict.items():
        options.append({
            "label": key,
            "value": key,
        })
    
    return html.Div([
        html.H6("Change the value in the text box to see callbacks in action!"),
        dcc.Checklist(
            id = "data-checklist",
            options = options
        ),
        html.Br(),
        html.Div(id="data-selected")
    ])

def create_gene_search():
    # column names should be the same as the ones from the dataframe containing the gene
    df = pd.DataFrame({"genes": [], "dataset": []})
    return  html.Div([
                html.H6("Find a gene in the dataset"),
                html.Div(["Gene name: ",
                        dcc.Input(id='gene-input', value="", type='text')]),
                html.Br(),
                dash_table.DataTable(
                    id='found-dataset',
                    columns=[{"name": i, "id": i} for i in df.columns],
                    data=df.to_dict('records'),
                    style_cell={'textAlign': 'center','min-width':'50px', 'width':'200px'})
            ])

def metadata_menu(metadata_df):
    return html.Div(children=[
        html.H6("Select sample options"),
        dcc.Dropdown(
            id="sample-dropdown",
            options=[{"label": i, "value": i} for i in metadata_df["Sample"]],
            value="Select samples",
            multi=True
        )
        # dropdown for tissue
        # dropdown for differentiation
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

    @dash_app.callback(
        Output(component_id='found-dataset', component_property="data"),
        [Input(component_id='gene-input', component_property='value')])
    def search_gene(user_input):
        genes_found, datasets = "", ""
        found_in_df = pd.DataFrame()
        if (user_input.strip()): 
            found_in_df = combined_genes[combined_genes["genes"].str.contains(user_input, case = False, na=False)]
            genes_found = " ".join(found_in_df["genes"].astype(str))
            datasets = " ".join(found_in_df["dataset"].astype(str))
        # return 'Gene(s) found: {} in dataset(s): {}'.format(genes_found, datasets)
        return found_in_df.to_dict('records')


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
        create_dfs_tick_box(data_dict),
        metadata_menu(data_dict["00_JBU_sequencing_metadata.tsv"]["df"]),
        create_gene_search() 
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






