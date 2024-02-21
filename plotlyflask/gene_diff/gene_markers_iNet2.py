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
        path = '/Users/vlad/Documents/Code/York/iNet/NB/Network_v2/healthy_hSBM/Top_50_ModCon_Rank.tsv'
    else: 
        path = '/Users/vlad/Documents/Code/York/iNet/NB/Network_v2/healthy_hSBM/Top_50_median.tsv'
        label = 'median'

    markers_df = pd.read_csv(path, sep='\t', dtype_backend='pyarrow')
    for col in markers_df:
        custom_traces.append({'genes': list(markers_df[col].values), 'title': f'{col}_{label}'})
  
    return custom_traces

