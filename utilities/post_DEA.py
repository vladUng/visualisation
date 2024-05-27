#!/usr/bin/env python
# -*-coding:utf-8 -*-
'''
@File    :   post_DEA.py
@Time    :   2023/01/20 14:08:16
@Author  :   Vlad Ungureanu
@Version :   1.0
@Contact :   vlad.ungureanu@york.ac.uk
@Desc    :   Script that process the output from Sleuth and prepares for JBU's visualisation tool
'''
import pandas as pd 
import numpy as np

from os import walk, path, makedirs

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

def prep_for_volcano(tcga_tpm_df, base_path, results_path, info_file, output_file, save_file = False):
    """
    Function that creates the file necessary for the volcano and scatter plots in the Visualisation tool.
    
    It receives a DataFrame with all the TPM values, which is used to create the DataFrame for the Viz tool by the following process:
    1. Gets the TPM values and applied log2(TPM+1)
    2. From the resultant df we keep only the genes used in the sleuth analysis (hence the need to read the sleuth output file)
    3. Calculate the average and median TPM values
    4. Compute the fold change based on the average and seperately fro median
    5. Add the data from sleuth to the output file

    Note: The median fold change is used for scatter plot while the average fold change is used for the volcano plot.

    Args:
        tcga_tpm_df (DataFrame): The TPM data, un-processed and NO log2 transformed. 
        mapping_cols (dictionary): Dictionary used to remap the samples by removing the suffix of -01A/-01B
        base_path (string): Where the file is saved 
        results_path (string): _description_
        info_file (string): _description_
        output_file (string): _description_
        save_file (bool, optional): _description_. Defaults to False.

    Returns:
        DataFrame: Resulting DataFrame
    """
    sleuth_results = pd.read_csv(results_path, sep="\t", names=["Genes", "p-value", "q-value"])
    pd_for_diff = pd.read_csv(base_path + info_file, sep="\t")

    # Create the new DataFrame and apply log2(TPM+1)
    dummy_df = pd.concat([pd.DataFrame(tcga_tpm_df["genes"]), pd.DataFrame(np.log2(tcga_tpm_df.iloc[:, 1:] + 1))], axis=1)

    # Add the missing values 
    sleuth_results = pd.concat([sleuth_results.set_index('Genes'), dummy_df.set_index('genes')], axis=1).fillna(1)
    sleuth_results['q-value'] = sleuth_results['q-value'].fillna(1)
    sleuth_results['p-value'] = sleuth_results['p-value'].fillna(1)
    # Bellow is dropping the genes that are found only in DEA and not in the expressed genes
    sleuth_results.dropna(inplace=True)

    # just keep the sample columns
    df = sleuth_results[dummy_df.columns[1:]].transpose().copy(deep=True)
    df.index.names = ["sample"]

    # set the cluster
    df["cluster"] = pd_for_diff.set_index("sample")["express"]

    print("Are arrays in sync? {}".format(np.array_equal(df.iloc[:, -1].reset_index(), pd_for_diff[["sample", "express"]])) )

    # calculate the median and avg TPM values
    fold_change = pd.DataFrame(df.columns[:-1], columns=["genes"])
    cluster_labels = pd_for_diff["express"].unique()
    new_labels = []
    for label in cluster_labels:
        new_labels.append("cluster_{}".format(label))
        fold_change[new_labels[-1] + "_med"] = df[df["cluster"] == label].iloc[:, :-1].median().values
        fold_change[new_labels[-1]] = df[df["cluster"] == label].iloc[:, :-1].mean().values

    # compute the fold change
    fold_change["fold_change_med"] = fold_change.iloc[:, 1] - fold_change.iloc[:, 3]
    fold_change["fold_change"] = fold_change.iloc[:, 2] - fold_change.iloc[:, 4]


    # assign the cluster labels
    fold_change["group"] = new_labels[0]
    fold_change.loc[fold_change["fold_change"] < 0, "group"] = new_labels[1]
    fold_change["-log10(q)"] = -np.log10(sleuth_results["q-value"])

    fold_change.set_index("genes", inplace=True)
    
    # Add the data from sleuth to the output file
    fold_change["q"] = sleuth_results["q-value"]
    fold_change["p"] = sleuth_results["p-value"]
    fold_change["pi"] = fold_change["-log10(q)"] * fold_change["fold_change"]

    fold_change.reset_index(inplace=True)
    fold_change.rename(columns={"index":"genes"}, inplace=True)
    
    if save_file:
        fold_change.to_csv(base_path + output_file,  index=False, sep="\t")
        
    return fold_change


# Inputs
version = "v1"

base_path = "../data/sel_prun/v1/"
viking_output = path.join(base_path, "viking/")
cluster_label = "dendrogram_cut"

# Read the data
tpm_df = pd.read_csv("../data/sel_prun/tum_TPMs_selected_genes_gc42_all_v4.tsv", sep="\t").rename(columns={'gene':'genes'})

raw_files = next(walk(viking_output), (None, None, []))[2]
experiments = [file.split("_results")[0] for file in raw_files]
# experiments.remove(".DS_Store")

dfs = {}
master_df = pd.DataFrame()
sel_cols = ["genes", "group", "pi", 'fold_change', '-log10(q)', "exp"]
for exp in experiments:
    results_path = f"{viking_output}/{exp}_results.tsv"
    info_file = f"{base_path}/info_files/{exp}.info"
    output_file = f"{viking_output}/{exp}_vulcano_labels.tsv"
    
    df = prep_for_volcano(tpm_df, '', results_path, info_file, output_file, save_file=True)
    df["exp"] = exp
    master_df = pd.concat([master_df, df[sel_cols]], axis=0)
    dfs["_".join(exp.split("_")[:-1])] = df
