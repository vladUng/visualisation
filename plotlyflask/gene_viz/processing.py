import time

import numpy as np
import pandas as pd

from plotlyflask.gene_viz import shared


def update_all_datasets(base_path):
    """
    Function which takes the names of all the files from the master.tsv files and then updates the all tsv file

    Args:
        base_path ([String]): Path to /data folder

    Returns:
        [type]: [description]
    """
    start_time = time.time()
    # open the master .csv file, containing the contents of the folder
    master_cv = pd.read_csv(base_path + "master.tsv", delimiter="\t")

    # iterate through the csv file and create a list of dictionaries with the following: <key> - name of the file, <value> object which contains a dataframe and description from the .csv file.
    all_tsv = pd.DataFrame()
    ret_dict = {}
    for _, row in master_cv.iterrows():
        filename, _ = row["Filename"], row["Description"]
        df = pd.read_csv(base_path + filename, delimiter="\t")
        print("Dataset: {} Shape {}".format(filename, df.shape))

        if "metadata" not in filename:
            if not all_tsv.empty:
                all_tsv = pd.merge(all_tsv, df, how="outer")
            else:
                all_tsv = df

    print("Finished loading the data in {}".format(time.time() - start_time))
    ret_dict["all_tsv"] = all_tsv

    start_time = time.time()
    all_tsv.to_csv(base_path + "all_data.tsv", sep="\t", index=False)
    print("Update .tsv file with the new data in {}".format(time.time() - start_time))


def process_data(df_metadata, all_tsv):
    """Process the metadata and the all_tsv file. The are two main operations happening:
    1. Replacing in specified cols the nan values and "?" with Not-known
    2. Fakes the _incl columns with values of "Y" to be as datasets. This involved to duplicate both the samples in metadata and in the `all_data.tsv`. We've marked the samples that are duplicated by the following rule "_incl_" + "including column initials".
    Args:
        df_metadata ([DataFrame]): The medata
        all_tsv ([DataFrame]): TPM dataframe
    Returns:
        (DataFrame, DataFrame): The two dataframes updated
    """

    df_metadata = df_metadata.replace("?", np.nan)
    # cols to replace nan with NA string
    cols = ["NHU_differentiation", "Tissue", "Gender", "Diagnosis", "Substrate"]
    for col in cols:
        df_metadata[col] = df_metadata[col].fillna("NA")

    # add the incl columns to the datasets and duplicate the samples
    incl_cols = [col for col in df_metadata.columns if "incl_" in col]
    new_rows_metadata = pd.DataFrame()
    new_cols_all_tsv = pd.DataFrame()
    for incl_col in incl_cols:
        # select samples
        metadata_selected = shared.filter_data(df_metadata, incl_col, ["Y"])
        selected_samples = metadata_selected["Sample"].values
        # get initials of the cols
        initials = "_incl_"
        for initial in incl_col.split("_"):
            initials += initial[0]
        # update the metadata
        duplicated_samples = [sample + initials for sample in metadata_selected["Sample"].values]
        metadata_selected = metadata_selected.drop("Sample", axis=1)
        metadata_selected["Sample"] = duplicated_samples
        metadata_selected["Dataset"] = incl_col
        metadata_selected["subset_name"] = incl_col
        new_rows_metadata = pd.concat([new_rows_metadata, metadata_selected], axis=0)

        # update all_tsv
        all_tsv_selected = all_tsv[selected_samples]
        all_tsv_selected.columns = duplicated_samples
        new_cols_all_tsv = pd.concat([new_cols_all_tsv, all_tsv_selected], axis=1)

    df_metadata = pd.merge(df_metadata, new_rows_metadata, how="outer")
    df_metadata["isTER"] = True
    for subtype in df_metadata["subset_name"].unique():
        if df_metadata[df_metadata["subset_name"] == subtype]["TER"].isnull().all():
            df_metadata.loc[df_metadata["subset_name"] == subtype, "isTER"] = False

    all_tsv = pd.concat([all_tsv, new_cols_all_tsv], axis=1)
    return df_metadata, all_tsv


##### Utilities functions #####


def correlations(df, gene_A, gene_B, isLog=False):
    spearman_df = df.corr(method="spearman").round(6)
    pearson_df = df.corr(method="pearson").round(6)

    pearson_data, spearman_data = [], []

    for idx, row in pearson_df.iterrows():
        if isLog:
            pearson_data.append({"corr_metric": "Log10({} + 1)".format(idx), "gene_a": row[gene_A], "gene_b": row[gene_B]})
        else:
            pearson_data.append({"corr_metric": idx, "gene_a": row[gene_A], "gene_b": row[gene_B]})

    for idx, row in spearman_df.iterrows():
        if isLog:
            spearman_data.append({"corr_metric": "Log10({} + 1)".format(idx), "gene_a": row[gene_A], "gene_b": row[gene_B]})
        else:
            spearman_data.append({"corr_metric": idx, "gene_a": row[gene_A], "gene_b": row[gene_B]})

    return pearson_data, spearman_data
