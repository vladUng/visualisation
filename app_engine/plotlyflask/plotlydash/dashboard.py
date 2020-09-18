import numpy as np
import pandas as pd
import dash
import dash_table
import dash_html_components as html
import dash_core_components as dcc
from plotlyflask.plotlydash.data import create_dataframe
from plotlyflask.plotlydash.layout import html_layout

import plotly.express as px

# filepaths and filenames
base_path = "/Users/vlad/Documents/Code/York/BU_clustering/src/data/"

bu_raw = "BU_TPMs.tsv"
# processed
bu_seq = 'BU_TPMs_sequencer-BC.tsv'
#metadata
bu_metadata = "200319_sequenced_sample_metadata.tsv"
# what is machine is used for sequencing are used
bu_machine = "BU_sequencer_batch.tsv"

# Differentiate tissue type after batch
all_diff_bc_raw = "all_diff_bc.tsv"
some_diff_bc_raw = "some_diff_bc.tsv"
po_bc_raw = "p0_bc_2.tsv"


def init_dashboard(server):
    """Create a Plotly Dash dashboard."""
    dash_app = dash.Dash(
        server=server,
        routes_pathname_prefix='/dashapp/',
        external_stylesheets=[
            'static/dist/css/styles.css',
            'plotlyflask/static/dist/css/styles.css',
            'https://fonts.googleapis.com/css?family=Lato'
        ]
    )

    # Load DataFrame
    df = create_dataframe()

    # Custom HTML layout
    dash_app.index_string = html_layout

    # Create Layout
    dash_app.layout = create_test_plot()
    return dash_app.server

def create_data_table(df):
    """Create Dash datatable from Pandas DataFrame."""
    table = dash_table.DataTable(
        id='database-table',
        columns=[{"name": i, "id": i} for i in df.columns],
        data=df.to_dict('records'),
        sort_action="native",
        sort_mode='native',
        page_size=300
    )
    return table

def create_test_plot():

    import plotlyflask.viz_tools.main_viz as viz

    fig1, fig2 = viz.test_bu_plots()

    return html.Div( children=[
           dcc.Graph(id="example-graph", figure=fig1),
           dcc.Graph(id="example-graph-2", figure=fig2)
        ],
        id='dash-container'
    )

    # df = pd.DataFrame({
    # "Fruit": ["Apples", "Oranges", "Bananas", "Apples", "Oranges", "Bananas"],
    # "Amount": [4, 1, 2, 2, 4, 5],
    # "City": ["SF", "SF", "SF", "Montreal", "Montreal", "Montreal"]
    # })

    # fig = px.bar(df, x="Fruit", y="Amount", color="City", barmode="group")
    # return dcc.Graph(
    #     id='example-graph',
    #     figure=fig
    # )

