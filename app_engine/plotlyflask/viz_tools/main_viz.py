import numpy as np
import pandas as pd 
import plotlyflask.viz_tools.utilities as utilities
import plotly.express as px


#loading the files, filepaths and filenames
base_path = "/Users/vlad/Documents/Code/York/BU_clustering/src/data/"

bu_raw = "BU_TPMs.tsv"
bu_seq = 'BU_TPMs_sequencer-BC.tsv' # processed
bu_metadata = "200319_sequenced_sample_metadata.tsv" # metadata
bu_machine = "BU_sequencer_batch.tsv" # the machine is used for sequencing are used

# Differentiate tissue type after batch
all_diff_bc_raw = "all_diff_bc.tsv"
some_diff_bc_raw = "some_diff_bc.tsv" #without Y2319 which hasn't differentiated properly (see BU_cluster_notebook)
po_bc_raw = "p0_bc_2.tsv" #only p0

# dataframes objects
df_bu_seq = pd.read_csv(base_path + bu_seq, sep='\t')
df_metadata = pd.read_csv(base_path + bu_metadata, sep='\t')
df_bu_machine = pd.read_csv(base_path + bu_machine, sep='\t')
df_bu_raw = pd.read_csv(base_path + bu_raw, sep='\t')
all_diff_bu_raw = pd.read_csv(base_path + all_diff_bc_raw, sep='\t')
some_diff_bu_raw = pd.read_csv(base_path + some_diff_bc_raw, sep='\t')
p0_bu_raw = pd.read_csv(base_path + po_bc_raw, sep='\t')

# 99 colours to choose from
colours_pool = px.colors.qualitative.Bold + px.colors.qualitative.D3 + px.colors.qualitative.G10 + px.colors.qualitative.T10 + px.colors.qualitative.Pastel + px.colors.qualitative.Prism + px.colors.qualitative.Vivid + px.colors.qualitative.Alphabet


#1
diff_bu_prcsd = utilities.filter_data(some_diff_bu_raw, th = 1.5, no_expressed = 3)
selected_samples = diff_bu_prcsd.columns.values[1:]

p0_bu_prcsd = utilities.filter_data(p0_bu_raw, th = 1.5, no_expressed = 3)
p0_samples = p0_bu_raw.columns.values[1:]

# 2
diff_set, p0_set = set(diff_bu_prcsd["genes"]), set(p0_bu_prcsd["genes"])
diff_surplus_gene =  diff_set - p0_set
p0_surplus_genes = p0_set - diff_set
common_genes = list(diff_set & p0_set)

diff_bu_prcsd = diff_bu_prcsd[diff_bu_prcsd["genes"].isin(common_genes)]
diff_bu_prcsd.reset_index(drop=True, inplace=True)
p0_bu_prcsd = p0_bu_prcsd[p0_bu_prcsd["genes"].isin(common_genes)]
p0_bu_prcsd.reset_index(drop=True, inplace=True)

# 3
diff_samples = utilities.ret_diff_samples([], selected_samples) #no samples to remove

#4
p0_sample_set = set([(p0_sample).split("-")[0] for p0_sample in p0_samples])
diff_sample_set = set([(diff_sample).split("-")[0] for diff_sample in diff_samples])

common_samples = list(p0_sample_set & diff_sample_set)
common_samples_p0 = [common_sample+"-P0" for common_sample in common_samples]
common_samples_diff = [common_sample+"-D" for common_sample in common_samples]

# 1. compute stats relative to the genes
var, median, std = diff_bu_prcsd.var(axis=1),  diff_bu_prcsd.median(axis=1),  diff_bu_prcsd.std(axis=1)
diff_bu_prcsd["variance"], diff_bu_prcsd["median"], diff_bu_prcsd["std"] = var, median, std
diff_bu_prcsd["r_median_std"] = diff_bu_prcsd["median"] /diff_bu_prcsd["std"]

# 2 - prepare the data 
diff_bu_prcsd.sort_values(by=["variance"] , ascending=False, inplace=True)
diff_bu_prcsd_2 = diff_bu_prcsd
diff_bu_prcsd_2["r_median_std"] = round(diff_bu_prcsd_2["median"] /diff_bu_prcsd_2["std"])
r_median_std_df = diff_bu_prcsd_2.sort_values(by=["r_median_std", "variance"] , ascending=False)


def test_bu_plots():
    # experiment configuration
    exp_config = {"name": "Cluster size 3","n_clusters": 3}
    config = (diff_samples, df_metadata, colours_pool)
    percent = 0.52

    hover_data = ["samples", "diagnosis", "ter_avg", "ter_sd", "labels_cystometric_capacity", "labels_date",
                "labels_ter_avg", "labels_ter_sd"]


    output_df, metrics_result = utilities.selected_genes_experiment(r_median_std_df, percent, config,
                                                    exp_config = exp_config, verbose=0, pca=True,
                                                    show_elbow = True, do_metrics=True, exp_id="_medv")
    fig1 = px.scatter(output_df, x = "tsne_2d_one", y= "tsne_2d_two", color = "Birch", text = "diagnosis", 
                    size="ter_avg", hover_data = hover_data, 
                    title = "The genes with the highest Median/variance ratio - Birch with t-SNE " + str(percent))
    fig2 = px.scatter(output_df, x = "tsne_2d_one", y= "tsne_2d_two", color = "RawKMeans", 
                    text = "samples", size="ter_avg",  hover_data =  hover_data,  
                    title = "The genes with the highest Median/variance ratio - Birch with t-SNE 3 clusters")

    
    return fig1, fig2