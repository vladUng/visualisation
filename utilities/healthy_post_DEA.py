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

# Pre-processing TPMs
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

def select_genes(tpm_df, no_genes = 3347, relative_selection = True):
    """
     It selects the most relative varied genes in the given DataFrame for the given number

    Args:
        tcga_tpm_df ([Dataframe]): The dataframe from where to select
        no_genes_selected (int, optional): [Genes to select]. Defaults to 3347.

    Returns:
        [type]: [description]
    """
    # dummy_df = pd.concat([pd.DataFrame(tcga_tpm_df["genes"]), pd.DataFrame(np.log2(tcga_tpm_df.iloc[:, 1:] + 1))], axis=1)
    dummy_df = np.log2(tpm_df.set_index("genes") + 1).reset_index()

    # remove all the genes w/ that have a lower expression value from `th` in `>10%` across the samples
    good_th = 0.05
    dummy_df = filter_data(dummy_df, th=np.log2(1.5), at_least_good_cols=dummy_df.shape[1]*good_th, idx_cols=["genes"])

    # make all the  values float
    dummy_df.set_index("genes", inplace=True)

    # acros samples
    print("####### Gene selection, num genes: {} #######".format(no_genes))
    if relative_selection:
        print("The genes selected by the highest standard deviation/median ration.")
        dummy_df["std"] = dummy_df.std(axis=1) / dummy_df.median(axis=1)
    else:
        print("The genes selected by the highest standard deviation; approached used by Robertson et al.")
        dummy_df["std"] = dummy_df.std(axis=1)

    most_varied_genes = list(dummy_df.sort_values(by="std", ascending=False).iloc[:no_genes].index)
    return most_varied_genes

def create_map_cols(tpm_df):
    """
     Remove the -01B and -01A - this needs to be run only once

    Args:
        tcga_tpm_df ([DataFrame]): where to remove

    Returns:
        [Dict]: Dictionary of the old vs new col name
    """
    mapping_cols = {}
    mapping_cols["genes"] = "genes"
    for col in tpm_df.columns.values[1:]:
        mapping_cols[col] = "-".join(col.split("-")[:-1])
    return mapping_cols

# Labeling functions

def encode_sleuth(row):
    if row["group"] == "cluster_0":
        return "LumP"
    elif row["group"] == "cluster_1":
        return "LumInf_NS"
    elif row["group"] == "cluster_2":
            return "Large_BaSq"
    elif row["group"] == "cluster_3":
        return "Small_BaSq"
    elif row["group"] == "cluster_4":
        return "NE-like"
    else:
        return "NaN"

def decode_sleuth(row):
    if row["comp_with"] == "LumP":
        return "c_0"
    elif row["comp_with"] == "LumInf_NS":
        return "c_1"
    elif row["comp_with"] == "Large_BaSq":
            return "c_2"
    elif row["comp_with"] == "Small_BaSq":
        return "c_3"
    elif row["comp_with"] == "NE-like":
        return "cluster_4"
    else:
        return "NaN"

def encode_sleuth_numeric(row):
    return int(row["group"].split("_")[1])

def get_most_sig(row):
    labels = row["exp"].split("_vs_")
    labels.remove(row["cluster"])
    return labels[0]

def apply_order_sig(df):
    # Goes through the name of the files(exp), splits the names in the two cluster labels (e.g ["Mixed"] '_vs_' ["Small_Ba_Sq_v4"]) and keeps ony the second (the first one we already know). Than we do some pre-processing
    for gene in df["genes"].unique()[:]:
        # order_sig = []
        order_sig = "-".join(df[df["genes"]==gene]["comp_with"].values)
        df.loc[df["genes"]==gene, "order_sig"] = order_sig

    return df 

def prep_for_volcano(tcga_tpm_df, mapping_cols, base_path, results_path, info_file, output_file, save_file = False, cluster_label='express'):
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

    # Difference between sleuth
    print(f"Diff between sleuth and all GE: {len(set(sleuth_results['Genes']) - set(tcga_tpm_df['genes']))}")
    print(f"Diff between GE and sleuth: {len(set(tcga_tpm_df['genes']) - set(sleuth_results['Genes']))}")

    # Select only the genes used in Sleuth
    df = dummy_df[dummy_df["genes"].isin(sleuth_results["Genes"])]
    # df.rename(columns=mapping_cols, inplace=True)
    df = df[["genes"] +  list(pd_for_diff["sample"].values)]
    df = df.set_index("genes").transpose()
    df.index.names = ["sample"]

    # set the cluster
    df["cluster"] = pd_for_diff.set_index("sample")[cluster_label]

    print("Are arrays in sync? {}".format(np.array_equal(df.iloc[:, -1].reset_index(), pd_for_diff[["sample", cluster_label]])) )

    # calculate the median and avg TPM values
    fold_change = pd.DataFrame(df.columns[:-1], columns=["genes"])
    cluster_labels = pd_for_diff[cluster_label].unique()
    new_labels = []
    for label in cluster_labels:
        new_labels.append(f"{label}")
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
    sleuth_results.set_index("Genes", inplace=True)
    
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

base_path = "/Users/vlad/Documents/Code/York/iNet/NB/Network_v2/test_iNet/Viking/"
viking_output = path.join(base_path, f"{version}")
# viz_tool = path.join(base_path, "Diff_exp/Viking/viz_tool/")
# cluster_label = "RawKMeans_CS_5"

# Read the data
tpm_df = pd.read_csv(f"{base_path}/{version}/healthy_data_all_gc42_v4.tsv", sep="\t").rename(columns={'gene':'genes'})


raw_files = next(walk(f"{viking_output}/results/"), (None, None, []))[2]
experiments = [file.split("_results")[0] for file in raw_files]
experiments
# experiments.remove(".DS_Store")

# experiments = ["Small_BaSq_vs_Lum_Inf_NS_v5"] 

dfs = {}
master_df = pd.DataFrame()
sel_cols = ["genes", "group", "pi", 'fold_change', '-log10(q)', "exp"]
for exp in experiments:
    if exp == '.DS_Store':
        continue
    print(f"\n###### {exp} ######")
    results_path = f"{viking_output}/results/{exp}_results.tsv"
    info_file = f"{viking_output}/info/{exp}.info"
    output_file = f"{viking_output}/results/labels/{exp}_v4_vulcano_labels.tsv"
    
    df = prep_for_volcano(tpm_df, {}, "", results_path, info_file, output_file, save_file=True, cluster_label='tissue_type')
    df["exp"] = exp
    # df = df.loc[df["genes"].isin(most_varied_genes)]
    master_df = pd.concat([master_df, df[sel_cols]], axis=0)
    dfs["_".join(exp.split("_")[:-1])] = df
