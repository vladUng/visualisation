#!/usr/bin/env python
# -*-coding:utf-8 -*-
'''
@File    :   prep_DEA.py
@Time    :   2023/01/20 09:16:00
@Author  :   Vlad Ungureanu
@Version :   1.0
@Contact :   vlad.ungureanu@york.ac.uk
@Desc    :   Takes the location of a .tsv which contains a DataFrame of samples with their cluster labelling. It outputs the .info file required for Viking
'''

import pandas as pd 
from itertools import combinations
from collections import Counter
from os import path, makedirs 

def prepare_for_viking(base_path, file_name, tpm_df, cluster="RawKMeans_CS_5"):
    """
    Function the pre-process a DataFrame for running Differentially Expressed Analysis with sleuth.

    Args:
        base_path (string): Where do you want to have the output file to be saved
        file_name (string): the name of the file
        tpm_df (DataFrame): The TPM dataframe
        cluster (str, optional): The column name of the cluster label. Defaults to "RawKMeans_CS_5".

    Returns:
        DataFrame: The resulting DataFrame.
    """
    df = pd.DataFrame(tpm_df[["Sample", cluster]].values, columns=["sample", "express"])

    df["path"] = "./01_h5_files/" + df["sample"] + "-01A_abundance.h5"
    # In TCGA there are some patients where they have been sampled or their sample was sequenced twice (1st-A, 2nd-B). We want to keep the latest (B)
    df.loc[df["sample"] == "TCGA-BL-A0C8", "path"] = "./01_h5_files/TCGA-BL-A0C8-01B_abundance.h5"
    df.loc[df["sample"] == "TCGA-BL-A13I", "path"] = "./01_h5_files/TCGA-BL-A13I-01B_abundance.h5"
    df.loc[df["sample"] == "TCGA-BL-A13J", "path"] = "./01_h5_files/TCGA-BL-A13J-01B_abundance.h5"
    df.loc[df["sample"] == "TCGA-GV-A3QK", "path"] = "./01_h5_files/TCGA-GV-A3QK-01B_abundance.h5"
    df.loc[df["sample"] == "TCGA-K4-A3WU", "path"] = "./01_h5_files/TCGA-K4-A3WU-01B_abundance.h5"
    df.loc[df["sample"] == "TCGA-K4-A4AB", "path"] = "./01_h5_files/TCGA-K4-A4AB-01B_abundance.h5"

    df.to_csv(base_path + file_name, index=False, sep="\t")
    return df

# Define the inputs
### base_path - input; base_path_results - output
base_path = "../../results/Stage I/gc42/"
filename = "VU_clustering.tsv"
cluster_label = "RawKMeans_42_CS_5"
version = "v5.1"

# Mapping; Sleuth only accepts integers 
labels_value = {0:"LumP", 1:"LumInf_NS", 2:"Large_BaSq",  3:"Small_BaSq", 4: "Mixed"}

# create the apropiate subfolders
base_path_results = path.join(base_path, "Diff_exp/Viking/" + version + "/")
if not path.exists(base_path_results):
    makedirs(base_path_results)

# Read the data
outputs = pd.read_csv(base_path + filename, sep="\t")

counter_values = Counter(outputs[cluster_label]) #for verifying
unique_values = outputs[cluster_label].unique()

unique_values.sort()

for comb in list(combinations(unique_values, 2)):
    cluster_1, cluster_2 = labels_value[comb[0]], labels_value[comb[1]]
    filename = "{}_vs_{}_v5.info".format(cluster_1, cluster_2)
    print("###Combinations {} for {}. Test the number of samples for each group:".format(comb, filename))

    to_save = pd.concat([outputs.loc[outputs[cluster_label] == comb[0]], outputs.loc[outputs[cluster_label] == comb[1]]])
    
    # We just check that the number of samples is the same in both cases 
    to_check = Counter(to_save[cluster_label])
    total_before = counter_values[comb[0]] +  counter_values[comb[1]]
    total_after = to_check[comb[0]] + to_check[comb[1]]
    if total_before != total_after and not to_save.isnull().values.any():
        print("❌ {}".format(filename))
        print("NaN values in the DataFrame", to_save.isnull().values.any())
        print("Before: Samples for {} are {} and for {} are {}".format(comb[0], counter_values[comb[0]], comb[1], counter_values[comb[1]]))
        print("After: Samples for {} are {} and for {} are {}".format(comb[0], to_check[comb[0]], comb[1], to_check[comb[1]]))
    else:
        print("✅ {}".format(filename))

    df = prepare_for_viking(base_path_results, file_name=filename, tpm_df=to_save.reset_index(), cluster=cluster_label)