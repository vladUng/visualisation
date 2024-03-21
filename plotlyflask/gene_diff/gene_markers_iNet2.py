#!/usr/bin/env python
# -*-coding:utf-8 -*-
'''
@File    :   gene_markers_iNet2.py
@Time    :   2024/02/21 08:55:40
@Author  :   Vlad Ungureanu
@Version :   1.0
@Contact :   vlad.ungureanu@york.ac.uk
@Desc    :   This file contains the gene markers for communities found with hSBM thorugh iNet2.
'''
import pandas as pd

def add_com_markers_raw(custom_traces, data_type='ModCon_Rank'):
    label = 'ModCon'
    if data_type == 'ModCon_Rank': 
        path = '/Users/vlad/Documents/Code/York/iNet/NB/Network_v2/healthy_hSBM_v4/Top_50_ModCon_Rank_v4.tsv'
    else: 
        path = '/Users/vlad/Documents/Code/York/iNet/NB/Network_v2/healthy_hSBM_v4/Top_50_median_v4.tsv'
        label = 'median'

    markers_df = pd.read_csv(path, sep='\t', dtype_backend='pyarrow')
    markers_df.columns= markers_df.columns.astype(int)
    markers_df = markers_df[markers_df.columns.sort_values()]
    sel_coms = [24, 26, 28, 30, 25, 22, 1, 29, 19]
    # sel_coms = all_coms
    for col in markers_df.columns:
        # if col not in sel_coms:
        #     continue
        custom_traces.append({'genes': list(markers_df[col].values), 'title': f'{col}_{label}'})
  
    return custom_traces

