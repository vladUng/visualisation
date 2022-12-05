#!/usr/bin/env python
# -*-coding:utf-8 -*-
'''
@File    :   pi.py
@Time    :   2022/11/11 13:14:23
@Author  :   Vlad Ungureanu
@Version :   1.0
@Contact :   vlad.ungureanu@york.ac.uk
@Desc    :   Creates the pi plot shown in the DEA page.

A pi plot has the following axis:
    * X: -log10(q) * fold_change <-- For comparison 1
    * Y: -log10(q) * fold_change <-- For comparison 2

The two comparisons are the results from the sleuth's DEA (see the Readme.md) file for the format of such file. 
'''


import numpy as np
import pandas as pd
import plotly.express as px

from plotlyflask.gene_diff import shared


def draw_pi_plot(df_1, df_2, filename_1, filename_2, selected_data, selected_genes):
    """
        Creates the pi plot

    Args:
        df_1 (DataFrame): The dataframe for the 1st comparison
        df_2 (DataFrame): Then dataframe for the 2nd comparison
        filename_1 (string): Name to be displayed on the plot of the 1st comparison
        filename_2 (string): Name to be displayed on the plot of the 2nd comparison
        selected_data (dict): List of the points selected
        selected_genes (_type_): List of the genes selected

    Returns:
        dict: The Pi figure
    """

    title = ""
    # Decide which is bigger, the number of genes may differ
    if df_1.shape[0] < df_2.shape[0]:
        first_df = df_2.copy(deep=True).set_index("genes")
        second_df = df_1.copy(deep=True).set_index("genes")
        title = "X: {} on  vs Y: {}".format(filename_2, filename_1)
    else:
        first_df = df_1.copy(deep=True).set_index("genes")
        second_df = df_2.copy(deep=True).set_index("genes")
        title = "X: {} on  vs Y: {}".format(filename_1, filename_2)
        
    # compute the values
    first_df["x"] = -np.log10(first_df["q"]) * first_df["fold_change"]
    second_df["y"] = -np.log10(second_df["q"]) * second_df["fold_change"]

    # Note the genes number may differ and we're setting the index of the DataFrame which has the most genes (i.e. rows)
    #  However, there might be some genes in the second df which are not in the first one. Regardless, we set the nan values to 0 (points will be added to the center)
    first_df.rename(columns={"group": "comp_1"}, inplace=True)
    second_df.rename(columns={"group": "comp_2"}, inplace=True)

    dummy_df = pd.concat([first_df[["x", "comp_1"]], second_df[["y", "comp_2"]]], axis=1).fillna(0).reset_index().rename(columns={"index":"genes"})
    
    dummy_df["main_colour"] = "PI_plot"
    # show in the hover everything apart from the main_colour column
    fig = px.scatter(dummy_df, x="x", y="y", hover_data=dummy_df.columns, color="main_colour", title=title)

    fig.update_traces(marker=dict(size=10, opacity=0.4), selector=dict(mode='markers'))
    
    x_col = "x"
    y_col = "y"

    # Add the selected points
    if selected_data:
        ranges = selected_data['range']
        selection_bounds = {'x0': ranges['x'][0], 'x1': ranges['x'][1],
                            'y0': ranges['y'][0], 'y1': ranges['y'][1]}

        # There is something weird going on with indexes
        selected_genes = [gene["customdata"][0] if "customdata" in gene.keys() else gene["text"] for gene in selected_data["points"]]

        selected_idxs = dummy_df[dummy_df["genes"].isin(selected_genes)].index

        fig.update_traces(selectedpoints=selected_idxs,
            mode='markers+text', 
            unselected={'marker': { 'opacity': 0.3 }
            })
    
        fig.add_shape(dict({'type': 'rect',
                        'line': { 'width': 1, 'dash': 'dot', 'color': 'darkgrey' } },
                       **selection_bounds))
    else:
        offset = 10
        selection_bounds = {'x0': np.min(dummy_df[x_col] - offset), 'x1': np.max(dummy_df[x_col] + offset),'y0': np.min(dummy_df[y_col] - offset), 'y1': np.max(dummy_df[y_col]) + offset}

    fig = add_anottations(first_df, second_df, dummy_df, fig)
    fig = show_selected_genes_pi(first_df.reset_index(), second_df.reset_index(), fig, selected_genes)

    return fig


def add_anottations(first_df, second_df, dummy_df, fig):
    offset = 10
    fig.add_shape(
        type='line',
        x0=dummy_df["x"].min()*1.5, y0=0,
        x1=dummy_df["x"].max()*1.5, y1=0,
        line=dict(color='Black', width=1),
        xref='x', yref='y'
    )

    fig.add_shape(
        type='line',
        x0=0, y0=dummy_df["y"].min()*1.5,
        x1=0, y1=dummy_df["y"].max()*1.5,
        line=dict(color='Black', width=1),
        xref='x',
        yref='y'
    )

    fig.add_annotation(
        showarrow=True,
        arrowhead=1,
        align = 'right',
        x=first_df["x"].max() + offset, y=0,
        text=first_df.loc[first_df["x"] == first_df["x"].max()]["comp_1"].values[0],
        opacity=0.7
    )

    fig.add_annotation(
        showarrow=True,
        arrowhead=1,
        align = 'left',
        x=first_df["x"].min() - offset, y= 0,
        text=first_df.loc[first_df["x"] == first_df["x"].min()]["comp_1"].values[0],
        opacity=0.7
    )
        
    fig.add_annotation(
        showarrow=True,
        arrowhead=1,
        align = 'right',
        y=second_df["y"].max() + offset, x=0,
        text=second_df.loc[second_df["y"] == second_df["y"].max()]["comp_2"].values[0],
        opacity=0.7
    )

    fig.add_annotation(
        showarrow=True,
        arrowhead=1,
        align = 'right',
        y=second_df["y"].min() - offset, x =0,
        text=second_df.loc[second_df["y"] == second_df["y"].min()]["comp_2"].values[0],
        opacity=0.7
    )

    return fig


def show_selected_genes_pi(df_1, df_2, fig, selected_genes):
    custom_traces = shared.create_custom_traces(selected_genes)
    colors =  px.colors.qualitative.Bold + px.colors.qualitative.Vivid 
    for idx, trace in enumerate(custom_traces): 
        fig.add_trace(shared.create_gene_trace(df_1, trace["genes"], name=trace["title"], marker_color=colors[idx], df_2=df_2))

    return fig 