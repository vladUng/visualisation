#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 14 15:07:04 2020

@author: vlad
"""
import numpy as np
import pandas as pd

################ Utilities functions ################

def filter_data(df, th, at_least_good_cols=3, idx_cols=["genes"]):
    """ Filter a given DataFrame by getting removing of the rows that have elements <= threshold.

    Args:
        df ([DataFrame]): The DataFrame from where we have to filter the data. This is has the genes as columns and samples as rows
        th ([float]): Threshold value of the unexpressed genes
        at_least_good_cols (int, optional): [The number of samples that fulfill the conditions]. Defaults to 3.

    Returns:
        [DataFrame]: DataFrame
    """
    # eliminating the first column
    df_prcsd = df.drop(idx_cols, axis=1)
    # compute the selected genes
    selected_genes_idxs = df_prcsd[df_prcsd >= th].dropna(thresh=at_least_good_cols).index.values
    selected_genes = df_prcsd.iloc[selected_genes_idxs]
    # add the genes names back
    cols = [df.loc[selected_genes_idxs, idx_col] for idx_col in idx_cols]
    cols.append(selected_genes)
    selected_genes = pd.concat(cols, axis=1)

    # reset indexes
    selected_genes.reset_index(drop=True, inplace=True)
    return selected_genes

def select_genes(tcga_tpm_df, no_genes = 3347, relative_selection = True):
    """
     It selects the most relative varied genes in the given DataFrame for the given number

    Args:
        tcga_tpm_df ([Dataframe]): The dataframe from where to select
        no_genes_selected (int, optional): [Genes to select]. Defaults to 3347.

    Returns:
        [type]: [description]
    """
    # dummy_df = pd.concat([pd.DataFrame(tcga_tpm_df["genes"]), pd.DataFrame(np.log2(tcga_tpm_df.iloc[:, 1:] + 1))], axis=1)
    dummy_df = np.log2(tcga_tpm_df.set_index("genes") + 1).reset_index()

    # remove all the genes w/ that have a lower expression value from `th` in `>10%` across the samples
    dummy_df = filter_data(dummy_df, th=np.log2(1.5), at_least_good_cols=dummy_df.shape[1]*0.9, idx_cols=["genes"])

    # make all the  values float
    dummy_df.set_index("genes", inplace=True)

    # acros samples
    print("####### Gene selection, num genes: {} #######".format(no_genes))
    if relative_selection:
        print("The genes selected by the highest standard deviation/median ration. So we choose the genes with the highest relative variance.")
        dummy_df["std"] = dummy_df.std(axis=1) / dummy_df.median(axis=1)
    else:
        print("The genes selected by the highest standard deviation; approached used by Robertson et al.")
        dummy_df["std"] = dummy_df.std(axis=1)

    most_varied_genes = list(dummy_df.sort_values(by="std", ascending=False).iloc[:no_genes].index)
    return most_varied_genes

def create_map_cols(tcga_tpm_df):
    """
     Remove the -01B and -01A - this needs to be run only once

    Args:
        tcga_tpm_df ([DataFrame]): where to remove

    Returns:
        [Dict]: Dictionary of the old vs new col name
    """
    mapping_cols = {}
    mapping_cols["genes"] = "genes"
    for col in tcga_tpm_df.columns.values[1:]:
        mapping_cols[col] = "-".join(col.split("-")[:-1])
    return mapping_cols
