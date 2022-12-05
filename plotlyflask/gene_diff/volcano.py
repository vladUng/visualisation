#!/usr/bin/env python
# -*-coding:utf-8 -*-
'''
@File    :   volcano.py
@Time    :   2022/11/11 13:48:23
@Author  :   Vlad Ungureanu
@Version :   1.0
@Contact :   vlad.ungureanu@york.ac.uk
@Desc    :   Creates the volcano plot in the DEA page

The plot is created from the sleuth's resultant file (see Readme for more info about it.)

A volcano plot has the following axis:
    * X: log2(fold_change) <-- Shows the change in the gene expression between teh two group
    * Y: log2(q) <-- Shows the significance of the change

'''


import dash_bio as dashbio
import numpy as np
import plotly.express as px

from plotlyflask.gene_diff import shared


def draw_volcano(df, fold_changes, selected_data, selected_genes):
    """ 
    Drawing the volcano plot

    # Todo: Check the below statements
     - There is something weird going on with indexes because we're setting customData to df['genes'].
     - It can get confusing as the data in df['genes'] will be put in alphabetical order and not correlated with the [x,y] points
     - Thus the real value of the genes will be different and is passed through selected_genes

    Args:
        df (dataframe): _description_
        fold_changes (dataframe): _description_
        selected_data (list): _description_
        selected_genes (list): the list of selected genes

    Returns:
        fig: figure object from plotly
    """
    fig = dashbio.VolcanoPlot(
        dataframe=df,
        effect_size="fold_change",
        gene="genes",
        snp=None,
        p="q",
        genomewideline_value=2.5,
        logp=True,
        effect_size_line=fold_changes,
        ylabel="-log10(q)",
        xlabel="log2(FC)",
        col='#2A3F5F',
        point_size=10,
        effect_size_line_width=4,
        genomewideline_width=2,
        highlight=True, 
        annotation="group"
    )
    
    x_col = "fold_change"
    y_col = "-log10(q)"

    # show the selected points
    if selected_data:
        ranges = selected_data['range']
        selection_bounds = {'x0': ranges['x'][0], 'x1': ranges['x'][1],
                            'y0': ranges['y'][0], 'y1': ranges['y'][1]}

        # finding the selected point
        selected_points = [gene["customdata"] for gene in selected_data["points"]]
        selected_idxs = df[df["genes"].isin(selected_points)].index
        fig.update_traces(selectedpoints=selected_idxs, mode='markers', unselected={'marker': { 'opacity': 0.3 } })
    else:
        selection_bounds = {'x0': np.min(df[x_col] - 2), 'x1': np.max(df[x_col] + 2),
                            'y0': np.min(df[y_col] - 2), 'y1': np.max(df[y_col]) + 2}
        fig.update_traces(marker=dict(size=10, opacity=0.4), selector=dict(mode='markers'))


    fig = show_selected_genes_vulcano(df, selected_genes, fig)

    fig.update_traces(customdata=df["genes"], textposition="top right")
    fig.update_layout(dragmode='select')

    fig.add_shape(dict({'type': 'rect', 'line': { 'width': 1, 'dash': 'dot', 'color': 'darkgrey' } }, **selection_bounds))

    return fig

def show_selected_genes_vulcano(df, selected_genes, fig):
    custom_traces = shared.create_custom_traces(selected_genes=selected_genes)
    colors =  px.colors.qualitative.Bold + px.colors.qualitative.Vivid 
    
    for idx, trace in enumerate(custom_traces): 
        
        fig.add_trace(shared.create_gene_trace(df, trace["genes"], name=trace["title"], marker_color=colors[idx]))

    return fig
