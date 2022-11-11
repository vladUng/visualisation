#!/usr/bin/env python
# -*-coding:utf-8 -*-
'''
@File    :   scatter.py
@Time    :   2022/11/11 13:37:54
@Author  :   Vlad Ungureanu
@Version :   1.0
@Contact :   vlad.ungureanu@york.ac.uk
@Desc    :   File creating the scatter plot as a complementary figure to the Volcano plot.

The plot is created from the sleuth's resultant file (see Readme for more info about it.)

A scatter plot has the following axis:
    * X: log2(fold_change) <-- Comparison 1
    * Y: log2(fold_change) <-- Comparison 2

'''


import numpy as np 

from plotlyflask.gene_diff import shared 
import plotly.express as px
import plotly.graph_objects as go

def plt_scatter(df, selected_points, known_markers=False):
    """
    Creates teh scatter plot

    Args:
        df (DataFrame): Thee data to be plotted
        selected_points (string): The points selected from the volcano plot
        known_markers (bool, optional): To display or not the custom traces. Defaults to False.

    Returns:
        dict: The Scatter figure
    """

    # select the columns with fold change based on the median
    clusters = [col for col in df.columns if ("med" in col) and ("fold" not in col)]
    fig = px.scatter(df, x=clusters[0], y=clusters[1], hover_data=df.columns, color="significance",  color_discrete_sequence=['#636EFA', "grey", '#EF553B'])

    fig.update_traces(marker=dict(size=10, opacity=0.4), selector=dict(mode='markers'))
    
    fig = plt_fc_lines(fig, df, fc_values = [0, 1, 2, 4])

    if not selected_points.empty:
        fig.add_trace(go.Scatter(x=selected_points[clusters[0]], y=selected_points[clusters[1]], mode="markers+text", text=selected_points["genes"], hoverinfo='all', textposition="top right", name="Selected Points"))

    if known_markers:
        custom_traces = shared.create_custom_traces()
        colors =  px.colors.qualitative.Vivid + px.colors.qualitative.Bold
        marker_size = 12
        for idx, trace in enumerate(custom_traces): 
            selected_df = df[df["genes"].isin(trace["genes"])]

            markers = {"size": marker_size, "color": selected_df.shape[0] * [colors[idx]], "symbol": "x"}
            trace = dict(type='scatter', x=selected_df[clusters[0]], y=selected_df[clusters[1]],  showlegend=True, marker=markers, text=selected_df["genes"], mode="markers+text" , name=trace["title"],  textposition="top right")

            fig.add_trace(trace)

    fig.update_layout(
        xaxis = dict(
            tickmode = 'array',
            tickvals = np.log2([1, 10, 100, 1000, 10000]),
            ticktext = ["0", "1", "10", "100", "1000", "10000"]
        ),
        yaxis = dict(
            tickmode = 'array',
            tickvals = np.log2([1, 10, 100, 1000, 10000]),
            ticktext = ["0", "1", "10", "100", "1000", "10000"]
        )
    )
    
    return fig

def draw_scatter(df, selected_genes):
    """
    Initiates all the function calls to draw the scatter plot

    Args:
        df (DataFrame): The data to be plotted
        selected_genes (list): Genes selected on volcano plots

    Returns:
        dict: The scatter figure
    """

    df = add_sig(df)
    selected_points = df[df["genes"].isin(selected_genes)]

    fig = plt_scatter(df, selected_points, True)
    fig.update_layout(clickmode='event+select')

    return fig

### Others
def plt_fc_lines(fig, df, fc_values):
    """

    On the given figure and dataframe, it shows the lines for the given fold change values

    Args:
        fig (dict): Plotly figure
        df (DataFrame): Data used on the figure
        fc_values (list): The list of the fold changes to display the fold change lines

    Returns:
        dict: The figure with the resultant lines
    """

    # get the maximum values across fold_change and the cluster groups
    med_cols = [col for col in df.columns if "med" in col and "fold" not in col]
    groups = ["fold_change"] + list(med_cols)
    max_value = round(df[groups].max().max())

    colours = ["Black", "Green", "goldenrod", "Red", "orange", "Green", "Purple"]
    line_type = [None, None, "dot", "dash", "dash"]
    for value in fc_values:

        if value:
            color = colours[value]
        else:
            color = "Black"

        line_dict = dict(color=color, width=3,  dash=line_type[value])

        fig.add_shape(type='line',
                x0=value, y0=0,
                x1=max_value, y1=max_value-value,
                line=line_dict,
                xref='x',yref='y', name="log2(FC)={}".format(value))

        if value:
            color = colours[value]
        else:
            color = "Black"

        line_dict = dict(color=color, width=3,  dash=line_type[value])

        fig.add_shape(type='line',
                x0=0, y0=value,
                x1=max_value-value+1, y1=max_value+1,
                line=line_dict,
                xref='x',yref='y', name="log2(FC)={}".format(value))

    return fig

# Classify the genes with significant changes
def significant_genes (row, labels):
   if row['fold_change_med'] < 1.0 and row['fold_change_med'] > -1.0:
      return 'Non-significant'
   elif row["fold_change_med"] > 1.0:
       return "Significant {}".format(labels[0])
   else: 
        return "Significant {}".format(labels[1])

def add_sig(df):
    """ 
    Functions that classifies the genes if are significant to one of the groups comparede in the DEA.

    Args:
        df (DataFrame): The data used

    Returns:
        DataFrame: The dataframe with the extra-col w/ indificating the type of significance
    """
    labels = df.columns.values[-2:]
    df["significance"] = df.apply(lambda row: significant_genes(row, labels), axis=1)
    return df