import pandas as pd
import numpy as np


from plotlyflask.viz_tools import utilities

base_path_results =  "/Users/vlad/Documents/Code/York/BU_clustering/results/Stage I/"
base_path_data = "/Users/vlad/Documents/Code/York/BU_clustering/data/TCGA/"

info_path = "/Users/vlad/Documents/Code/York/visualisation/visualisation/data/ScatterPlot/infoFiles/"

def get_tpms_df():
    tcga_tpm_df = pd.read_csv(base_path_data + "Remap/TCGA-BLCA_gc38_filtered-remap_TPMs.tsv", sep="\t")
    tcga_metadata_df = pd.read_csv(base_path_data + "metadata_tcga.csv")

    # remove the samples that are duplicates (i.e. resamples and are marked by having 01A#).
    # the good ones are the 01B
    samples_to_remove = []
    for sample in tcga_tpm_df.columns[1:-1]:
        sample_type = sample.split("-")[3]
        if len(sample_type) > 3:
            samples_to_remove.append(sample)

    # drop from metadata, tpm
    tcga_tpm_df.drop(samples_to_remove, axis=1, inplace=True)
    tcga_metadata_df.drop(samples_to_remove, axis=1, inplace=True)

    selected_genes = utilities.select_genes(tcga_tpm_df, no_genes=3500)

    return tcga_tpm_df, selected_genes

def significant_genes (row, labels):
   if row['FC'] < 1 and row['FC'] > -1:
      return 'Non-significant'
   elif row["FC"] > 1:
       return "Significant {}".format(labels[0])
   else: 
        return "Significant {}".format(labels[1])

def prc_fc_comp(tcga_tpm_df, sleuth_results, exp, drop_outliers = False):

    # get the tpms
    pd_for_diff = pd.read_csv(info_path + exp + ".info", sep="\t")

    ##### Apply pre-processing by Robertson et al. #####
    dummy_df = pd.concat([pd.DataFrame(tcga_tpm_df["genes"]), pd.DataFrame(tcga_tpm_df.iloc[:, 1:])], axis=1)

    df = dummy_df[dummy_df["genes"].isin(sleuth_results["genes"])]

    mapping_cols = utilities.create_map_cols(tcga_tpm_df)
    df.rename(columns=mapping_cols, inplace=True)

    df = df[["genes"] +  list(pd_for_diff["sample"].values)]
    df = df.set_index("genes").transpose()
    df.index.names = ["sample"]
    df["cluster"] = pd_for_diff.set_index("sample")["express"]

    print("Are arrays in sync? {}".format(np.array_equal(df.iloc[:, -1].reset_index(), pd_for_diff[["sample", "express"]])) )

    fold_change = pd.DataFrame(df.columns[:-1], columns=["genes"])
    cluster_labels = pd_for_diff["express"].unique()
    new_labels = []
    for label in cluster_labels:
        new_labels.append("cluster_{}".format(label))
        fold_change[new_labels[-1]] = df[df["cluster"] == label].iloc[:, :-1].median().values

    fold_change[new_labels] = np.log2(fold_change[new_labels]+1)
    outliers = []
    if drop_outliers:
        outliers = fold_change.sort_values(by=new_labels, ascending=False).iloc[0:10]["genes"].values 
        fold_change = fold_change[~fold_change["genes"].isin(outliers)]

    fold_change["FC"] = fold_change[new_labels[0]] - fold_change[new_labels[1]]
    fold_change["sig"] = fold_change.apply(lambda row: significant_genes(row, new_labels), axis=1 )
    return fold_change, outliers
