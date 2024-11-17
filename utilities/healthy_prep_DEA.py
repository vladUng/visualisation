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

# def prepare_for_viking(base_path, file_name, tpm_df, cluster="RawKMeans_CS_5", cluster_label='dendrogram_label', h5_path="./01_h5_files/"):
#     """
#     Function the pre-process a DataFrame for running Differentially Expressed Analysis with sleuth.

#     Args:
#         base_path (string): Where do you want to have the output file to be saved
#         file_name (string): the name of the file
#         tpm_df (DataFrame): The TPM dataframe
#         cluster (str, optional): The column name of the cluster label. Defaults to "RawKMeans_CS_5".

#     Returns:
#         DataFrame: The resulting DataFrame.
#     """
#     df = pd.DataFrame(tpm_df[["Sample", cluster, cluster_label]].values, columns=["sample", "express", cluster_label])

#     df["path"] = h5_path + df["sample"] + "-01A/"+ "abundance.h5"
#     # In TCGA there are some patients where they have been sampled or their sample was sequenced twice (1st-A, 2nd-B). We want to keep the latest (B)
#     df.loc[df["sample"] == "TCGA-BL-A0C8", "path"] = f"{h5_path}/TCGA-BL-A0C8-01B/abundance.h5"
#     df.loc[df["sample"] == "TCGA-BL-A13I", "path"] = f"{h5_path}/TCGA-BL-A13I-01B/abundance.h5"
#     df.loc[df["sample"] == "TCGA-BL-A13J", "path"] = f"{h5_path}/TCGA-BL-A13J-01B/abundance.h5"
#     df.loc[df["sample"] == "TCGA-GV-A3QK", "path"] = f"{h5_path}/TCGA-GV-A3QK-01B/abundance.h5"
#     df.loc[df["sample"] == "TCGA-K4-A3WU", "path"] = f"{h5_path}/TCGA-K4-A3WU-01B/abundance.h5"
#     df.loc[df["sample"] == "TCGA-K4-A4AB", "path"] = f"{h5_path}/TCGA-K4-A4AB-01B/abundance.h5"

#     df['express'] = df['express'].astype(int)
#     # df.drop(columns=cluster_label, inplace=True)
#     df.to_csv(base_path + file_name, index=False, sep="\t")
#     return df

# Define the inputs
### base_path - input; 
# base_path_results - output
# base_path = "../../results/Stage I/gc42/"
base_path = '../data/cluster_analysis/'
filename = "VU_clustering_v4.tsv"
cluster_model = "RawKMeans_labels" #it needs to be numeric
cluster_label = 'KMeans_labels_6'
version = "v1"
h5_path = '/mnt/scratch/projects/biol-cancerinf-2020/Raw-Data/TCGA/BLCA/01_RNAseq/kallisto-gencode-v42/'

# Mapping; Sleuth only accepts integers 
# labels_values = {13:"LumP", 12:"LumInf", 4:"Large_BaSq",  5:"Small_BaSq", 3: "Mes-like"} //for Net_I
labels_values = {22: "Med_IFNG", 20: "Low_IFNG", 3: "High_IFNG", 0: "Lump", 4: "NE", 1: "Lum_InfNs"}

# create the apropiate subfolders
base_path_results = path.join(base_path + version + "/info_files/")
if not path.exists(base_path_results):
    makedirs(base_path_results)

# Read the data
outputs = pd.read_csv(base_path + filename, index_col='Sample', sep="\t")

counter_values = Counter(outputs[cluster_model]) #for verifying
unique_values = outputs[cluster_model].unique()

unique_values.sort()

for comb in list(combinations(unique_values, 2)):
    cluster_1, cluster_2 = labels_values[comb[0]], labels_values[comb[1]]
    filename = f"{cluster_1}_vs_{cluster_2}_{version}.info"
    print("###Combinations {} for {}. Test the number of samples for each group:".format(comb, filename))

    to_save = pd.concat([outputs.loc[outputs[cluster_model] == comb[0]], outputs.loc[outputs[cluster_model] == comb[1]]])
    
    # We just check that the number of samples is the same in both cases 
    to_check = Counter(to_save[cluster_model])
    total_before = counter_values[comb[0]] +  counter_values[comb[1]]
    total_after = to_check[comb[0]] + to_check[comb[1]]
    if total_before != total_after and not to_save.isnull().values.any():
        print("❌ {}".format(filename))
        print("NaN values in the DataFrame", to_save.isnull().values.any())
        print("Before: Samples for {} are {} and for {} are {}".format(comb[0], counter_values[comb[0]], comb[1], counter_values[comb[1]]))
        print("After: Samples for {} are {} and for {} are {}".format(comb[0], to_check[comb[0]], comb[1], to_check[comb[1]]))
    else:
        print("✅ {}".format(filename))

    df = prepare_for_viking(base_path_results, file_name=filename, tpm_df=to_save.reset_index(), cluster=cluster_model, cluster_label = cluster_label,h5_path=h5_path)