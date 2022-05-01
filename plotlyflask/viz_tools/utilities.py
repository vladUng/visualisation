#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 14 15:07:04 2020

@author: vlad
"""
import numpy as np
import pandas as pd

# for plotting
import plotly.graph_objects as go
import plotly.express as px
from plotly.subplots import make_subplots
from scipy.cluster.hierarchy import dendrogram

from scipy.spatial.distance import pdist, squareform
from sklearn.preprocessing import StandardScaler, Normalizer
from sklearn.decomposition import PCA
from sklearn.preprocessing import LabelEncoder
from sklearn.manifold import TSNE

from sklearn import cluster, mixture
from sklearn.cluster import KMeans
from sklearn.neighbors import kneighbors_graph

from sklearn import metrics

import warnings
import itertools

from os import path
import pickle
import time
import random

from collections import Counter
# from yellowbrick.cluster import KElbowVisualizer


################ Utilities functions ################

def ret_diff_samples(samples_to_remove, all_samples):
    """Function that remove samples from the given array.

    Args:
        samples_to_remove ([String])
        all_samples ([String])

    Returns:
        [String]: the selected samples
    """
    selected_samples = []
    for elem in all_samples:
        if not elem in samples_to_remove:
            selected_samples.append(elem)
    return selected_samples

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


def add_prcsd_metdata(df, df_meta):
    """ Adds to a df metadata columns with information such as TER avg, std, tissue type etc.
        Note: It doesn't all the metadata, but some selected, which was needed.

    Args:
        df ([DataFrame]): Has the rows as samples and the columns either PC or t-sne values. Also, it's the one to whom the metadata is appended it
        df_meta ([DataFrame]): Df that has the samples metadat, like age, TER etc

    Returns:
        [LabelEncoder]: This it's returned as it may be useful to classify later.
        [list]: Column with all the metadata added
    """
    diff_seq, tissue_types, diagnosis = [], [], []
    ter_avgs, ter_sds= [], []
    labels_ter_avg, labels_ter_sd, labels_date, labels_cystometric = [], [], [], []

    #remove the white characters be aware this will change the argument passed
    df_meta["sample"] = df_meta["sample"].str.strip()
    # Add TER avg and std, but first we have to impute some values
    df_meta["TER_avg"].fillna(df_meta["TER_avg"].median(), inplace=True)
    df_meta["TER_sd"].fillna(df_meta["TER_sd"].median(),  inplace=True)
    # add labels
    label_encoder = LabelEncoder()
    df_meta["labels_date"] = label_encoder.fit_transform(df_meta['date'])
    df_meta["labels_cystometric_capacity"] = label_encoder.fit_transform(df_meta['cystometric_capacity'])

    # add the sequencing machine column, N/A if the tissue type doesn't fit
    for sample_name in df["samples"]:
        sample_id = sample_name.split("-")[0].strip()
        tissue_type = sample_name.split("-")[1].strip()
        # tissue types
        tissue_types.append(tissue_type)
        # diff_seq
        sample = df_meta[df_meta["sample"] == sample_id]
        if len(sample["diff_sequencer"].values) > 0:
            if tissue_type == "D": # seperate for diff and P0 tissue
                diff_seq.append(sample["diff_sequencer"].values[0])
            elif tissue_type == "P0":
                diff_seq.append(sample["P0_sequencer"].values[0])
            else:
                diff_seq.append("N/A")
        else:
            diff_seq.append("N/A")

        # ter avg and standard deviation
        ter_avg = sample["TER_avg"].values[0]
        ter_sd = sample["TER_sd"].values[0]
        # labels for ter_avg
        if ter_avg < 100:
            labels_ter_avg.append("absent")
        elif ter_avg > 100 and ter_avg < 500:
            labels_ter_avg.append("weak")
        elif ter_avg > 500:
            labels_ter_avg.append("tight")

        # labels for ter_sd
        if ter_sd < 100:
             labels_ter_sd.append("absent")
        elif ter_sd > 100 and ter_sd < 500:
             labels_ter_sd.append("weak")
        elif ter_sd > 500:
             labels_ter_sd.append("tight")

        # add to the arrays
        ter_avgs.append(ter_avg)
        ter_sds.append(ter_sd)
        labels_date.append(sample["labels_date"].values[0])
        labels_cystometric.append(sample["labels_cystometric_capacity"].values[0])
        # clinical diagnosis
        diagnosis.append(sample["clin_diagnosis"].values[0])

    # add the data to the given df
    df["tissue_type"] = tissue_types
    df["diff_sequencer"] = diff_seq
    df["diagnosis"] = diagnosis
    df["ter_avg"] = ter_avgs
    df["labels_ter_avg"] = labels_ter_avg
    df["ter_sd"] = ter_sds
    df["labels_ter_sd"] = labels_ter_sd
    df["labels_date"] = labels_date
    df["labels_cystometric_capacity"] = labels_cystometric

    # we may need this later to unpack some of the encoded values
    cols =  ["samples", "tissue_type",  "labels_ter_avg", "ter_avg", "labels_cystometric_capacity", "labels_date", "labels_ter_sd",  "ter_sd", "diff_sequencer"]
    return cols

def feature_scalling(df, selected_features, scaleType="standard_transpose"):
    """ Features scales the received DataFrame. There are the following scaling options implemented with the scikit:
        1. row-wise: l1, l2, max and standard_transpose
        2. col-wise: standard scaller

    Check scikit-learn docs for more information on how each works

    Args:
        df ([DataFrame]): The data to which to apply the operations, it contains the first column as the gene
        selected_features ([type]): We may not want to apply all the samples, this means that the returned df may have a different size
        scaleType (str, optional): Specify which type of scaler you want to apply. Defaults to "standard_transpose".

    Returns:
        [DataFrame]: The result from the scaling operation
    """
    dummy_df = df.loc[:, selected_features]
    if scaleType == "standard_transpose":
        scalled = StandardScaler().fit_transform(dummy_df.transpose().values)
    elif scaleType == "standard":
        scalled = StandardScaler().fit_transform(dummy_df.values)
    elif scaleType == "l2":
        scalled = Normalizer(norm='l2').fit_transform(dummy_df.values)
    elif scaleType == "l1":
        scalled = Normalizer(norm='l1').fit_transform(dummy_df.values)
    elif scaleType == "max":
        scalled = Normalizer(norm='max').fit_transform(dummy_df.values)
    return pd.DataFrame(scalled)

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

    # acros samples
    print("####### Gene selection, num genes: {} #######".format(no_genes))
    if relative_selection:
        print("The genes selected by the highest standard deviation/median ration. So we choose the genes with the highest relative variance.")
        dummy_df["std"] = dummy_df.std(axis=1) / dummy_df.median(axis=1)
    else:
        print("The genes selected by the highest standard deviation; approached used by Robertson et al.")
        dummy_df["std"] = dummy_df.std(axis=1)

    most_varied_genes = dummy_df.sort_values(by="std", ascending=False).iloc[:no_genes]["genes"]
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


################ PCA, elbow method, metrics, functions ################

def apply_pca(df, samples_names, genes_names, pca_components, transpose = True):
    """ For a given DataFrame we apply Principal Component Analysis

    Args:
        df (DataFrame): The data to be processed
        samples_names ([list]): Neeed to add back to the returned df
        genes_names ([pd.Series]): Neeed to add back to the returned df
        pca_components ([int]): number of components
        transpose (bool, optional): This is needed so that we know how to add the genes/samples back to the df. Defaults to True.

    Returns:
        [DataFrame]: The PCA with samples and genes
    """
    # apply PCA
    pca = PCA(n_components = pca_components)
    pca_bu = pca.fit_transform(df)
    print('Variation per principal component {}'.format(pca.explained_variance_ratio_))

    # generate columns labels
    pca_col_labels = []
    for idx in range(pca_components):
        pca_col_labels.append("PC_" + str(idx + 1)) # start from 1

    pca_bu_df = pd.DataFrame(data = pca_bu, columns = pca_col_labels)
    # add dataframes labels accordingly
    dummy_df = df.copy(deep=True)
    if transpose:
        dummy_df.columns = genes_names
        dummy_df.insert(0, "samples", samples_names)
        pca_bu_df.insert(0, "samples", samples_names)
    else:
        dummy_df.columns = samples_names
        dummy_df.insert(0, "genes", genes_names)
        pca_bu_df.insert(0, "genes", genes_names)

    return pca_bu_df, pca

def apply_kmeans(df_pca, centroids = 4):
    kmeans = KMeans(n_clusters = centroids)
    centroids_labels = kmeans.fit_predict(df_pca)
    return centroids_labels

def elbow_method(df, min_k, max_k):
    """ Apply elbow method on K-means for a given data (DataFrame)

    Args:
        df ([DataFrame]): Input data
        min_k ([Int]):
        max_k ([Int]):

    Returns:
        [plotly.Figure, list]: Return the figure to be plotted and the sum squared distances
    """
    sum_dist, K = [], range(min_k, max_k)
    for k in K:
        km = KMeans(n_clusters=k)
        km = km.fit(df)
        sum_dist.append(km.inertia_)

    data = go.Scatter(x=list(K), y=sum_dist, mode="lines+markers")
    layout = go.Layout( title= "Elbow method",
        xaxis=go.layout.XAxis(title="Clusters"),
        yaxis=go.layout.YAxis(title="Sum of squared distances"))
    fig = go.Figure(data=data, layout=layout)
    return fig, sum_dist

def scale_elbow_pca(df, selected_samples, no_pca, no_clusters, scaling_type):
    scaled_diff_bu = feature_scalling(df, selected_samples, scaling_type)
    if scaling_type == "standard_transpose":
        dummy_df = scaled_diff_bu
    else:
        dummy_df = scaled_diff_bu.transpose()

    pca_bu_df, _pca = apply_pca(df=dummy_df, samples_names=selected_samples, genes_names=df["genes"],pca_components=no_pca, transpose=True)
    fig, _ = elbow_method(pca_bu_df.iloc[:, 1:], 1, no_clusters)
    return fig

def compute_cluster_metrics(cluster_models, data, id=""):
    """ Calculates the ‘cityblock’, ‘cosine’, ‘euclidean’, ‘l1’, ‘l2’, ‘manhattan’ metrics for a given set of cluster models

    Args:
        cluster_models ([list]): List of cluster models such as Kmeans, Birch, Ward etc.
        data ([DataFrame]): The data used to trained the cluster models
        id ([String]): The id to append to the cluster names. This is useful for a set of experiments. Defaults to ""

    Returns:
        [DataFrame]: DataFrame used for storing the cluster metrics
    """

    clusters_metrics, cluster_names = [], []
    for name, cluster_model in cluster_models:
        if name != "GaussianMixture":
            labels = cluster_model.labels_
        else:
            labels = cluster_model.predict(data)

        # ‘cityblock’, ‘cosine’, ‘euclidean’, ‘l1’, ‘l2’, ‘manhattan’
        silhoute_metric = metrics.silhouette_score(data, labels, metric='euclidean')
        silhoute_metric_2 = metrics.silhouette_score(data, labels, metric='manhattan')
        silhoute_metric_3 = metrics.silhouette_score(data, labels, metric='cosine')
        calinski_harabasz_metric = metrics.calinski_harabasz_score(data, labels)
        davies_bouldin_metric = metrics.davies_bouldin_score(data, labels)

        clusters_metrics.append([silhoute_metric, silhoute_metric_2, silhoute_metric_3, calinski_harabasz_metric, davies_bouldin_metric])
        cluster_names.append(name + id)

    metrics_names = ["Silhoute_euclidean", "Silhoute_manhattan", "Silhoute_cosine", "Calinski_habrasz", "Davies_bouldin"]
    cluster_metrics = pd.DataFrame(np.around(clusters_metrics, decimals=5), columns = metrics_names)
    cluster_metrics.insert(0, "Cluster", cluster_names)

    return cluster_metrics

def compute_experiments(config, dfs, samples, colours_pool, verbose=0):
    """ This is an extension of the clustering methods function. This will plot the metrics for a received input_data from t_sne

    Args:
        config (touple): IT contains the configuration for the experiments: output_df, datasets, selected clusters and column names.
            - The output_df is the the data that is fed to the clustering methods function, which itself applies different clustering algorithms to the output_df
            - datasets represents an array of dictionary of the form: { (data), experiment_settings }. Where:
                - data is input_data,
                - experiment settings is the form of the default base in this function:
                    default_base = {'quantile': .2,
                                    'eps': .3,
                                    'damping': .5,
                                    'preference': 5,
                                    'n_neighbors': 5,
                                    'n_clusters': 3,
                                    'min_samples': 5,
                                    'xi': 0.05,
                                    'min_cluster_size': 0.1,
                                    "birch_th": 1.7,
                                    'name': ""
                                    }
            - selected_clusters are the the cluster models IDs to be used, currently are supported the following:
             "RawKMeans", "Ward", "Birch", "GaussianMixture", "SpectralClustering", "AgglomerativeClustering", "Gaussian Mixture", "Affinity Propagation", "DBSCAN", "OPTICS", "Hierarchical Clustering"

            - column names, axis labels for the plots, which has to be the columns of the input_data

        dfs ([touple]): Two DataFrames one containing the dataframe which will contain the last results (e.g. clustering labels of the experiments) and the other the metadata.

        samples ([list]): The list of samples to be used

        colours_pool ([list]): List of colours to be used for plots

        verbose (int, optional): If to show plots. Defaults to 0.

    Returns:
        [DataFrame, list]: The dataframe with the latests results from the experiments and all the results from the experiments
    """
    output_df, meta_df = dfs[0], dfs[1]
    # the default configurations
    default_base = {'quantile': .2,
                    'eps': .3,
                    'damping': .5,
                    'preference': 5,
                    'n_neighbors': 5,
                    'n_clusters': 3,
                    'min_samples': 5,
                    'xi': 0.05,
                    'min_cluster_size': 0.1,
                    "birch_th": 1.7,
                    'name': ""}

    # unfold the parameters
    input_data, datasets, selected_clusters = config[0], config[1], config[2]
    results = clustering_methods(datasets, default_base, samples, selected_clusters)
    if config[3] != None:
        col_names = config[3]

    # plotting the clusters
    for (name, output, cluster_models) in results:
        # added the aditional data
        algorithm_names = output.iloc[:, 1:].columns
        output_df["samples"] = samples
        output_df[algorithm_names] = output[algorithm_names]

        output_df["Birch"] = order_labels_size(output_df["Birch"])
        output_df["RawKMeans"] = order_labels_size(output_df["RawKMeans"])

        if verbose:
            # plotting clusters
            plot_clusters(output_df, algorithm_names, col_names[0], col_names[1], colours_pool, name, marker_text="diagnosis")
            # metrics and their plots
            metrics = compute_cluster_metrics(cluster_models, input_data)
            plot_cluster_metrics(metrics, algorithm_names, name)

    _ = add_prcsd_metdata(output_df, meta_df)
    return output_df, results

################ Plotting functions ################

def draw_heatmap(df, fullpath, name, save):
    """Creating a heatmap from a given DataFrame

    Args:
        df ([DataFrame]): the data to be used for heatmap
        fullpath ([String]): Where to save the figure
        name (): Name of the heatmap which will be used as filename
        save ([Bool]): to sve or not the figure
    """

    # get the x_labels which in our case are the genes name
    x_labels = df['genes']
    df.drop(['genes'], axis = 1, inplace=True)

    # get the y labels which are representing by the sample names
    y_labels = df.columns.tolist()
    data_dist = pdist(df[y_labels].values)

    data = [ go.Heatmap( z=squareform(data_dist),
            colorscale='spectral',
            x=y_labels,
            y=x_labels)]

    layout = go.Layout( title=name,
    xaxis=go.layout.XAxis(title="Samples"),
    yaxis=go.layout.YAxis(title="Genes"))

    fig = go.Figure(data=data, layout = layout)
    fig.show()
    if save:
        save_fig(fullpath + "/%s" % (name), fig)


def clustering_methods(datasets, default_base, samples, selected_clusters=None):
    """  This function is the core of applying different clustering methods to a dataset. It can run different experiments with different datasets and configuration for the algorithms. It's a modification of the scikit-learn [blogpost](https://scikit-learn.org/stable/modules/clustering.html). The following clusters are supported (name, id):

        1. Kmeans - RawKMeans"
        2. Mini Batch KMeans - MiniBatchKMeans
        3. Ward - Ward
        4. Birch - Birch
        5. Gaussian Mixture Models - GaussianMixture
        6. Affinity Propagation - AffinityPropagation
        7. SpectralClustering - Spectral Clustering
        8. DBSCAN - DBSCAN
        9. OPTICS - OPTICS
        10. Hierarchical Clustering - Hierarchical Clustering

    Args:
        datasets ([dic]): List of the touples which containes the datasets, which is a DataFrame and the cluster models parameters in the form of a dictionary that needs to be override. See default_base of the parameters that can be override.
        default_base ([dict]): The configurations to be override for an experiments. Defaults to 3 clusters and birch_th 1.7. In case it needs to be override, below are the acceptable parameters and defualt values:
                    {'quantile': .2,
                    'eps': .3,
                    'damping': .5,
                    'preference': 5,
                    'n_neighbors': 5,
                    'n_clusters': 3,
                    'min_samples': 5,
                    'xi': 0.05,
                    'min_cluster_size': 0.1,
                    "birch_th": 1.7,
                    'name': "Name of the experiment" }
        samples ([String]): The list of samples
        selected_clusters ([String], optional): List of the strings that are the cluster models supported. Defaults to None, which means that all the available cluster models will be used.

    Returns:
        [list]: List of touples which contains the following data (in  this order): the experiment name, cluster model output and the cluster model object itself.
    """

    # array which contains a list with the name of dataset, the cluster models and the output labels
    ret_val = []
    for (dataset, algo_params) in datasets:
        # update parameters with dataset-specific values
        params = default_base.copy()
        params.update(algo_params)

        X = dataset
        # connectivity matrix for structured Ward
        connectivity = kneighbors_graph(X, n_neighbors=params['n_neighbors'], include_self=False)
        # make connectivity symmetric
        connectivity = 0.5 * (connectivity + connectivity.T)

        # ============
        # Create cluster objects
        # ============
        two_means = cluster.MiniBatchKMeans(n_clusters=params['n_clusters'])
        kmeans = KMeans(n_clusters = params['n_clusters'], max_iter = 1000)
        ward = cluster.AgglomerativeClustering( n_clusters=params['n_clusters'], linkage='ward',  connectivity=connectivity)
        birch = cluster.Birch(n_clusters=params['n_clusters'], threshold=params["birch_th"])
        affinity_propagation = cluster.AffinityPropagation(damping=params['damping'], preference=params['preference'])
        gmm = mixture.GaussianMixture(n_components=params['n_clusters'], covariance_type='diag', max_iter = 500) #we know from prev experiments this is a good configuration
        spectral = cluster.SpectralClustering(n_clusters=params['n_clusters'], eigen_solver='arpack',affinity="nearest_neighbors")
        dbscan = cluster.DBSCAN(eps=params['eps'], min_samples = params['min_samples'])
        optics = cluster.OPTICS(min_samples=params['min_samples'], xi=params['xi'], min_cluster_size=params['min_cluster_size'])
        average_linkage = cluster.AgglomerativeClustering(linkage="average", affinity="cityblock", n_clusters=params['n_clusters'], connectivity=connectivity)
        hierarchical_clutering = cluster.AgglomerativeClustering(distance_threshold=0,  n_clusters=None,  linkage="average",  connectivity=connectivity)

        clustering_algorithms = (
            ("RawKMeans", kmeans),
            ('MiniBatchKMeans', two_means),
            ('Ward', ward),
            ('Birch', birch),
            ('GaussianMixture', gmm),
            ('AffinityPropagation', affinity_propagation),
            ('SpectralClustering', spectral),
            ('AgglomerativeClustering', average_linkage),
            ('DBSCAN', dbscan),
            ('OPTICS', optics),
            ('Hierarchical Clustering', hierarchical_clutering)
        )

        output_algorithms = pd.DataFrame(samples, columns=['samples'])
        output_models = []
        for name, algorithm in clustering_algorithms:
            # selected_clusters is None means that all are selected
            if selected_clusters is not None:
                if name not in selected_clusters:
                    continue #skip if it's not in the loop
            # catch warnings related to kneighbors_graph
            with warnings.catch_warnings():
                warnings.filterwarnings(
                    "ignore",
                    message="the number of connected components of the " +
                    "connectivity matrix is [0-9]{1,2}" +
                    " > 1. Completing it to avoid stopping the tree early.",
                    category=UserWarning)
                warnings.filterwarnings(
                    "ignore",
                    message="Graph is not fully connected, spectral embedding" +
                    " may not work as expected.",
                    category=UserWarning)
                algorithm.fit(X)

            if hasattr(algorithm, 'labels_'):
                y_pred = algorithm.labels_.astype(np.int)
            else:
                y_pred = algorithm.predict(X)

            output_algorithms[name] = y_pred
            output_models.append(( name, algorithm ))

        # this need refactoring to returning a dictionary with key, pairs as it will be clearer
        ret_val.append((params["name"], output_algorithms, output_models))

    return ret_val

def plot_pcs_two_three(df_pca, pca_bu, hover_data, dot_colour=None, dot_text=None, dot_size=None, fig_path = None):
    """ Plots two figures for a dataframe that is formed out of the Principal components + metadata. One of the plots is 2D for the 2 PC, the other is 3D for 3 components.

    Args:
        df_pca ([DataFrame]):  DataFrame that contains the 3 PCs and the metadata.
        pca_bu ([PCA]): The PCA object which holds information of the PCA such as expalined variance ratio
        hover_data ([list]): Data to be displayed when point is covered
        dot_colour ([String]): The coloring by a field from metadata. Defaults to None.
        dot_text ([String], optional): Any field from metadata. Defaults to None.
        dot_size ([String], optional): Any field from metadata. Defaults to None.
        fig_path (String, optional): If we want to save the figures we neeed to specify a path + name, to which we append the name of the field used for coloring the clusters. Defaults to None.
    """
    if dot_size == None:
        dot_size = "ter_avg"
    if dot_text == None:
        dot_text = "diagnosis"

    title_2pc = ("Colouring by %s" % dot_colour) + " 2PCs - PC_1({:.3f}%), PC_2({:.3f}%)".format(pca_bu.explained_variance_ratio_[0] * 100, pca_bu.explained_variance_ratio_[1]* 100)
    title_3pc = ("Colouring by %s" % dot_colour) + " 3 PCs -  PC_1({:.3f}%), PC_2({:.3f}%), PC_3({:.3f}%)".format( pca_bu.explained_variance_ratio_[0] * 100, pca_bu.explained_variance_ratio_[1] * 100, pca_bu.explained_variance_ratio_[2] * 100)
    fig1 = px.scatter(df_pca, x="PC_1", y="PC_2", color=dot_colour, text=dot_text, size=dot_size,
                      hover_data = hover_data, title=title_2pc)
    fig2 = px.scatter_3d(df_pca, x="PC_1", y="PC_2", z="PC_3", text=dot_text, size=dot_size,
                         color=dot_colour, hover_data = hover_data,
                         title= title_3pc)

    fig1.update_traces(textposition='bottom left')
    fig2.update_traces(textposition='bottom left')
    fig1.show()
    fig2.show()

    if fig_path != None:
        save_fig(dot_colour + "_2_PC", fig1, base_path=fig_path)
        save_fig(dot_colour + "_3_PC", fig2, base_path=fig_path)

def trace_2d(df, colours_pool, title, x_axis, y_axis, hover_cols=None, dot_colour=None, marker_text=None, marker_text_position=None):
    """ This creates a 2D trace that it can be later be used for subplots in plotly. Easy to customise, the only constriant is that the x and y axis has to be part of the df argument.

    Args:
        df ([DataFrame]): The data that contains the points to be plotted.
        colours_pool ([list]): The list of colours which are used for plotting
        title ([String]): The title of the plot
        x_axis ([String]): X label, it has to be part of the DataFrame
        y_axis ([String]): Y label, it has to be part of the DataFrame
        hover_cols ([list], optional): [description]. Defaults to (centroids and samples for backwards compatibility).
        dot_colour ([String], optional): The string which is the column name of the df and by that it's done the colouring. Defaults to centroids (for backwards compatibility).
        marker_text ([String], optional): The string which is a column name from the df. Defaults to None.
        marker_text_position ([String], optional): A string identified form plotly that's used to center accordingly the marker text. Defaults to None.

    Returns:
        [trace]: The trace object that was created
    """
    hover, colours, markers = [], [], {}
    mode_markers = "markers"
    text_list = []
    if hover_cols == None:
        hover_cols =  ["centroids", "samples", "diagnosis", "labels_ter_avg", "labels_ter_sd"]
    if dot_colour == None:
        dot_colour = "centroids"
    if marker_text != None:
        mode_markers = "markers+text"
        text_list = df[marker_text].values
        if marker_text_position == None:
            marker_text_position="bottom left"
    for _, row in df.iterrows():
        centroid = row[dot_colour]
        # create the hover data
        hover_string = ""
        for hover_col in hover_cols:
            hover_string += "<br>%s=%s" %(hover_col, str(row[hover_col]))
        hover_string +=  "<br>" + x_axis +"=%{x}" + "<br>" + y_axis +"=%{y}"
        hover.append(hover_string)

        colours.append(colours_pool[centroid])

    markers["color"] = colours
    markers["size"] = 6
    trace = dict(type='scatter', x=df[x_axis], y=df[y_axis],
           hovertemplate=hover, showlegend=True, name=title, marker=markers,
           mode=mode_markers, text=text_list, textposition=marker_text_position)
    return trace

def trace_3d(df, colours_pool, title, x_axis, y_axis, z_axis, hover_cols=None, dot_colour=None, marker_text=None, marker_text_position=None):
    """ This creates a 3D trace that it can be later be used for subplots in plotly. Easy to customise, the only constriant is that the x,y, z axis has to be part of the df argument.

    Args:
        df ([DataFrame]): The data that contains the points to be plotted.
        colours_pool ([list]): The list of colours which are used for plotting
        title ([String]): The title of the plot
        x_axis ([String]): X label, it has to be part of the DataFrame
        y_axis ([String]): Y label, it has to be part of the DataFrame
        z_axis ([String]): Z label, it has to be part of the DataFrame
        hover_cols ([list], optional): [description]. Defaults to (centroids and samples for backwards compatibility).
        dot_colour ([String], optional): The string which is the column name of the df and by that it's done the colouring. Defaults to centroids (for backwards compatibility).
        marker_text ([String], optional): The string which is a column name from the df. Defaults to None.
        marker_text_position ([String], optional): A string identified form plotly that's used to center accordingly the marker text. Defaults to None.


    Returns:
        [trace]: The trace object that was created
    """
    mode_markers = "markers"
    hover, colours, markers, text_list = [], [], {}, []
    if hover_cols == None:
        hover_cols = ["centroids", "samples", "diagnosis", "labels_ter_avg", "labels_ter_sd"]
    if dot_colour == None:
        dot_colour = "centroids"
    if marker_text != None:
        mode_markers = "markers+text"
        text_list = df[marker_text].values
        if marker_text_position == None:
            marker_text_position="bottom left"
    for _, row in df.iterrows():
        centroid = row[dot_colour]
        # create the hover data
        hover_string = ""
        for hover_col in hover_cols:
            hover_string += "<br>%s=%s" %(hover_col, str(row[hover_col]))
            hover_string +=  "<br>" + x_axis +"=%{x}" + "<br>" + y_axis +"=%{y}" + "<br>" + z_axis +"=%{z}"
            hover.append(hover_string)

            colours.append(colours_pool[centroid])

    markers["color"] = colours
    markers["size"] = 10
    trace = dict(type='scatter3d', x=df[x_axis], y=df[y_axis], z=df[z_axis],
           hovertemplate=hover,  showlegend=True, name=title, marker=markers,
           mode=mode_markers, text=text_list, textposition=marker_text_position)
    return trace


def plot_clusters(df, plot_cols, x_label, y_label, colours_pool, exp_name, hover_data=None, marker_text=None, marker_text_position = None):
    """ Function that plots the results from the clustering methods. It can received any number of columns values strings and it will display 3 on a row. This means that it will always display plots with 3 columns

    Args:
        df ([DataFrame]): At the moment only it supports 2d result DataFrames. This means can be either a t-sne or PCA
        plot_cols ([list]): A list of strings with the columns values to be ploted
        x_label([String]): The column of x which should be plotted, for t-sne it can be sne_2d_one. It has to be in df arg
        y_label ([String]): The column of y which should be plotted, for t-sne it can be sne_2d_two. It has to be in df arg
        colours_pool ([list]): The list of colours which are used for plotting
        exp_name ([String]): Name for experiment which is used for plots
        marker_text ([String], optional): The string which is a column name from the df. Defaults to None.
        marker_text_position ([String], optional): A string identified form plotly that's used to center accordingly the marker text. Defaults to None.
    """
    num_cols = 3
    num_rows = int(np.ceil(len(plot_cols)/num_cols))
    if hover_data == None:
        hover_data = [plot_cols, "samples"]

    fig = make_subplots(rows=num_rows, cols=num_cols, subplot_titles = plot_cols, horizontal_spacing=0.1, vertical_spacing=0.1)

    traces = []
    for plot_col in plot_cols:
        hover_data.append(plot_col)
        trace = trace_2d(df, colours_pool, plot_col, x_label, y_label, hover_cols=hover_data, dot_colour=plot_col, marker_text=marker_text, marker_text_position=marker_text_position)
        traces.append(trace)
        hover_data.remove(plot_col)

    # add the traces to the subplot
    subplot_titles = []
    idx_row, idx_col = 1,1
    for trace in traces:
        subplot_titles.append(trace["name"])
        fig.add_trace(trace, row=idx_row, col=idx_col)
        # What we do here, we increment the column and row idxs. This means we increment the column at each iteration and reset it when it is divisible with the max number of columnbs. After that we increment the row idx
        if idx_col % num_cols == 0:
            idx_col = 0
            idx_row +=1
        idx_col += 1

    layout = go.Layout(title_text=exp_name,
                        xaxis=go.layout.XAxis(title=x_label),
                        yaxis=go.layout.YAxis(title=y_label))
    fig.update_layout(layout, height=1000)
    fig.show()

def plot_cluster_metrics(metrics, algorithm_names, exp_name, hide_text = False):
    """ Creates the the figures for a given DataFrame of metrics.

    Args:
        metrics (DataFrame): Holds the metrics values, columns names being represented by the algorithm for which the metrics
        algorithm_names ([list]): List of strings of the algorithms that are being run
        exp_name ([String]): Name of the experiment
        hide_text ([Bool]): If True the points text is hidden. This is usefull when there are lots of points. Defaults to False
    """
    fig = go.Figure()
    # define x axis
    metrics_names = ["Silhoute_euclidean", "Silhoute_manhattan", "Silhoute_cosine", "Calinski_habrasz", "Davies_bouldin"]
    mode, text, traces = 'lines+markers+text', metrics["Cluster"].values, []
    random_x = np.linspace(0, 15, metrics.shape[0])
    if hide_text:
        mode="lines+markers"
    for metrics_name in metrics_names:
        trace = go.Scatter(x=random_x, y=metrics[metrics_name], mode= mode, text = text, hoverinfo='all', textposition="top right")
        traces.append(trace)

    fig = make_subplots(rows=3, cols=2,
                        subplot_titles = ["Silhouette euclidean (higher the better)",
                                          "Silhouette manhattan (higher the better)",
                                          "Silhouette cosine (higher the better)",
                                          "Calinski-Harabrasz (higher the better)",
                                          " Davies-Bouldin (lower the better)"],
                        shared_yaxes= False,
                        horizontal_spacing=0.1, vertical_spacing=0.1 )

    for idx, trace in enumerate(traces):
        fig.add_trace(trace, row=int(idx/2)+1, col=idx%2+1)

    layout = go.Layout( title_text="Different clustering metrics " + exp_name)
    fig.update_layout(layout)
    fig.show()

def genes_distrib(df, no_genes=100):
    """ Finds the distribution of the most abundant genes in a given dataframe and plots a histogram

    Args:
        df (DataFrame): The dataframe with the genes
        no_genes (int, optional): The number of the most abundant genes. Defaults to 100.
    """
    #get the gene names of the most abundant ones
    working_df = drop_columns_if_needed(df).copy(deep=True)
    working_df["variance"] = working_df.var(axis=1)
    working_df["r_median_std"] =  working_df.median(axis=1) / working_df.std(axis=1)
    working_df.sort_values(by="variance", ascending=False, inplace=True)
    most_abundant_genes = working_df.iloc[:no_genes, :]
    #find the indexes of these genes in the df, which is where we look for their distribution
    gene_idxs = (df[df["genes"].isin(most_abundant_genes["genes"].values)].index.values)
    most_abundant_genes["distrib"] = (gene_idxs * 100) / (working_df.shape[0])
    fig = px.histogram(most_abundant_genes, x="distrib", hover_data=["genes", "variance", "r_median_std",  most_abundant_genes.index], marginal="rug")
    fig.show()

################ % of the dataset experiments ################

def test_variance(data_df, data_ref, diff_samples, df_metadata, colours_pool, percent_range, exp_config=None):
    """ Function that runs a set of experiments with the goal to test how the most abundant of something (e.g. variance or std/ratio) is driving the clusters. This means that the received data (data_df) has to be sorted (descending, highest first) by the respective column (e.g. column). Using the given percent_range the function will iterated through them and eliminate gradually the lowest genes of that particular feature

    Args:
        data_df (DataFrame): Has to be sorted by the feature that it needs to be tested.
        data_ref (DataFrame): Reference for results.
        diff_samples ([list]): List of string samples that we're using in the DataFrame
        df_metadata ([DataFame]): The DataFrame that contains the metadata
        colours_pool ([list]): Colours pool for cluster colours
        percent_range ([range]): The range of percentages
        exp_config ([dict], optional):  exp_config ([dictionary], optional): The configurations to be override for an experiments. Defaults to 3 clusters and birch_th 1.7. In case it needs to be override, below are the acceptable parameters and defualt values:
                    {'quantile': .2,
                    'eps': .3,
                    'damping': .5,
                    'preference': 5,
                    'n_neighbors': 5,
                    'n_clusters': 3,
                    'min_samples': 5,
                    'xi': 0.05,
                    'min_cluster_size': 0.1,
                    "birch_th": 1.7,
                    'name': "Name of the experiment" }

    Returns:
        [DataFrame, DataFrame, dict]: The output of the clustersing in the experiments used, the changes in the cluster sizes and the cluster transition of the samples between experiments.
    """

    # needed to surpress some warnings about modifyin sliced Panda series
    pd.options.mode.chained_assignment = None
    # define the structures that wrere stored
    membership_transition, cluster_sizes = {}, pd.DataFrame()
    clusters_outputs = pd.DataFrame(data_ref[["samples", "Birch", "RawKMeans"]].values, columns = ["samples", "Birch-100_0", "RawKMeans-100_0"])
    clusters_outputs.insert(0, "change_count", np.zeros(data_ref["samples"].shape[0]))
    cluster_sizes["Birch-ref"], cluster_sizes["KMeans-ref"] = get_cluster_size(data_ref["Birch"]), get_cluster_size(data_ref["RawKMeans"])
    for percent in percent_range:
        percent = round(percent, 5) #sometimes arrange yields a large number decimals
        run_str = str(percent*100).replace(".","_")
        no_genes = int(round(data_df.shape[0] * percent))
        print("Out of the 14k {}% are selected = {}".format(percent*100, no_genes))

        # setup the data to be used
        diff_bu_sel = data_df.iloc[:no_genes, :]
        data_sel =  diff_bu_sel[diff_samples].transpose()

        # run exp
        datasets_sel, selected_clusters = [(data_sel, {"name": "Cluster size 3", "n_clusters": 3})], ["RawKMeans", "Birch"]
        if exp_config != None:
            datasets_sel = [(data_sel, exp_config)]
        config, dfs = (diff_bu_sel, datasets_sel, selected_clusters, None), (pd.DataFrame(), df_metadata)
        output_df_sel, _ = compute_experiments(config, dfs, diff_samples, colours_pool, verbose = 0)

        # get birch and kmeans from previous exp for comparison
        last_birch, _ = clusters_outputs.iloc[:, -2], clusters_outputs.iloc[:, -1]

        # track the changes of each individual samples from one cluster to other
        old_birch = last_birch[last_birch != output_df_sel["Birch"]]
        new_birch_change = output_df_sel[last_birch != output_df_sel["Birch"]][["samples", "Birch"]]
        birch_changes = pd.concat([new_birch_change, old_birch], axis = 1)

        birch_changes.columns = ["samples", "Birch" + run_str, last_birch.name]
        membership_transition["run-" + run_str] = birch_changes

        changes_count = np.where(output_df_sel["Birch"].values != clusters_outputs.iloc[:, -2].values)
        clusters_outputs["change_count"].loc[changes_count] += 1

        # add additional info for analysis
        cluster_sizes["Birch-"+run_str] = get_cluster_size(output_df_sel["Birch"])
        cluster_sizes["KMeans-"+run_str] = get_cluster_size(output_df_sel["RawKMeans"])
        clusters_outputs["Birch-"+run_str] = output_df_sel["Birch"]
        clusters_outputs["KMeans-"+run_str] = output_df_sel["RawKMeans"]

    # revert back
    pd.options.mode.chained_assignment = "warn"

    clusters_outputs = clusters_outputs.reindex(sorted(clusters_outputs.columns, reverse = True), axis=1)
    cluster_sizes = cluster_sizes.reindex(sorted(cluster_sizes.columns, reverse = True), axis=1)
    return clusters_outputs, cluster_sizes, membership_transition

def selected_genes_experiment(diff_bu, percent, config, exp_config=None, tsne_config=None, pca=True, show_elbow=False, do_metrics = False, exp_id = "", fig_path=None, verbose=0):
    """ Here we apply differet analysis methods to a given DataFrame. It starts by creating a t-sne object and then, depending on the parameters it starts applying:
        - PCA
        - clustering models
        - clustering metrics
        - elbow method & visualisation

    Args:
        diff_bu ([DataFrame]): Genes x Samples. this is important as we drop the genes column and transpose the df
        percent ([Double/Float]): The precentage of genes to be included, if you want to include all genes use 1.0
        config ([type]): A touple which has on [0] - samples, [1] - the metadata DataFrame, [2] - Colours pool for plots.
        exp_config ([dictionary], optional): The configurations to be override for an experiments. Defaults to 3 clusters and birch_th 1.7. In case it needs to be override, below are the acceptable parameters and defualt values:
                    {'quantile': .2,
                    'eps': .3,
                    'damping': .5,
                    'preference': 5,
                    'n_neighbors': 5,
                    'n_clusters': 3,
                    'min_samples': 5,
                    'xi': 0.05,
                    'min_cluster_size': 0.1,
                    "birch_th": 1.7,
                    'name': "Name of the experiment" }
        tsne_config ([dictionary], optional): A touple to override the t-sne settings. Default values are perplexity, n_iter, learning_rate = 20, 10000, 200. To change them pass a dictionary of the form:
                     { "perplexity":10,
                       "learning_rate": 15,
                       "n_iter": 10000 }
        pca (bool, optional): Applies PCA if True. Defaults to True.
        show_elbow (bool, optional): Showing the elbow methods . Defaults to False.
        do_metrics (bool, optional): If to calculate metrics. Defaults to False
        exp_id (String, optional): An identifier string that is appended to the clustering models used when we are computing the metrics. This is useful if we want to compare different metrics. Defaults to an empty stirng.
        fig_path ([String], optional): The path for figures. If nothing is specified we don't save the figures. Defaults to None.
        verbose (int, optional): Display all the plots & other useful information. Defaults to 0.

    Returns:
        [DataFrame, list]: The DataFrame of the t-sne with the metadata added and a list with the result of the cluster metrics
    """
    # unpack data
    samples, df_metadata, colours_pool = config[0], config[1], config[2]
    no_genes = int(round(diff_bu.shape[0] * percent))
    print("Out of the 14k {}% are selected = {}".format(percent*100, no_genes))

    # prepare the DataFrame that we are going to use
    diff_bu_sel = diff_bu.iloc[:no_genes, :]
    data_sel = drop_columns_if_needed(diff_bu_sel).drop(["genes"], axis=1).transpose()

    # build the t-sne object
    perplexity, n_iter, learning_rate = 20, 10000, 200
    if tsne_config != None:
        perplexity, n_iter, learning_rate = tsne_config["perplexity"], tsne_config["n_iter"], tsne_config["learning_rate"]
    col_names = ["tsne_2d_one", "tsne_2d_two"]
    tsne = TSNE(learning_rate=learning_rate, perplexity=perplexity, n_iter=n_iter, n_iter_without_progress=1000, init = 'pca')
    tsne_sel = tsne.fit_transform(data_sel)
    tsne_sel_df = pd.DataFrame(tsne_sel, columns=col_names)
    tsne_sel_df["samples"] = samples

    if pca:
        pca_bu_df, pca_bu = apply_pca(df= data_sel, samples_names=samples,genes_names=diff_bu_sel["genes"], pca_components=3, transpose=True)
        data_sel = pca_bu_df.iloc[:, 1:]
    else:
        print("###Runing non-PCA mode")

    # Define experiments
    datasets_dif = [(data_sel, {"name": "Cluster size 3","n_clusters": 3, 'birch_th': 1.7})]
    if exp_config != None:
        datasets_dif = [(data_sel, exp_config)]
    selected_clusters = ["RawKMeans", "Birch"]
    config = (data_sel, datasets_dif, selected_clusters, col_names)
    dfs = (tsne_sel_df, df_metadata)

    # run experiment and unpacked the results
    _ , results_pca = compute_experiments(config, dfs, samples, colours_pool, verbose = 0)
    name, output, cluster_models = results_pca[0][0], results_pca[0][1], results_pca[0][2]

    algorithm_names = output.iloc[:, 1:].columns
    tsne_sel_df["samples"] = samples
    tsne_sel_df[algorithm_names] = output[algorithm_names]
    # keep the colours relevant to their size
    tsne_sel_df["RawKMeans"] = order_labels_size(tsne_sel_df["RawKMeans"])
    tsne_sel_df["Birch"] = order_labels_size(tsne_sel_df["Birch"])
    _ = add_prcsd_metdata(tsne_sel_df, df_metadata)

    if pca:
        pca_bu_df[algorithm_names] = output[algorithm_names]
        pca_bu_df["RawKMeans"] = order_labels_size(pca_bu_df["RawKMeans"])
        pca_bu_df["Birch"] = order_labels_size(pca_bu_df["Birch"])
        cols = add_prcsd_metdata(pca_bu_df, df_metadata)
        if verbose:
            plot_pcs_two_three(pca_bu_df, pca_bu, cols, "RawKMeans", fig_path= fig_path)
            plot_pcs_two_three(pca_bu_df, pca_bu, cols, "Birch", fig_path= fig_path)
    # not very eficient to calculate this every time, but there are instances when the metrics results should be displayed differently
    metrics_results = []
    if do_metrics:
        metrics_results = compute_cluster_metrics(cluster_models, data_sel, "_%d_%s" % ((percent*100), exp_id))
        if verbose:
            plot_cluster_metrics(metrics_results, algorithm_names, name)

    if verbose and show_elbow:
        elbow_metrics = ["distortion", "silhouette", "calinski_harabasz"]
        for elbow_metric in elbow_metrics:
            km = KMeans()
            visualizer = KElbowVisualizer(km, k=(2,12), metric=elbow_metric)
            visualizer.fit(data_sel)
            visualizer.show()
    return tsne_sel_df, metrics_results

def drop_columns_if_needed(df):
    """Function which drops some extra columns that we have used to sort a DataFrame.

    Currently we drop the following: "mean", "variance","r_mean_std", "r_median_std", "median", "std"
    Args:
        df ([DataFrame]): The DataFrame from which we drop the columns

    Returns:
        [DataFrame]: The new DataFrame
    """
    to_drop = []
    extra_columns = ["mean", "variance","r_mean_std", "r_median_std", "median", "std" ]
    for extra_column in extra_columns:
        if extra_column in df.columns.values:
            to_drop.append(extra_column)

    return df.drop(to_drop, axis=1)

################ Random gene configurations experiments ################

def run_random_config(diff_bu, iterations, no_genes, samples):
    """ 
    Function that computes the random configuration experiments. This means that the received data (diff_bu) is put through a given number of iterations (iterations) to run birch and kmeans.

    Args:
        diff_bu (DataFrame): The dataset to be used for the clustering 
        iterations (Int): The number of iterations to run the birch & kmeans
        no_genes (Int): Number of integers
        samples (list): List of string samples that we're using in the DataFrame

    Returns:
        [Dictionary]: The dictionary with the results which includes the metrics for all iterations, the best 1000 configurations, best run etc.
    """

    import warnings
    warnings.filterwarnings("ignore")

    # start_time = time.time()
    diff_bu.reset_index(inplace=True,  drop=True)
    df_data = diff_bu[samples].transpose()
    df_data.reset_index(inplace=True, drop=True)

    best_genes = pd.DataFrame()
    best_run, best_silhouette, saved_config, scores, stats  = 0, 0, [], [], []
    start_time_2 = time.time()
    for k in range(0, iterations):
        random_genes = random.sample(range(0, df_data.shape[1]), no_genes)
        kmeans = KMeans(n_clusters = 3, max_iter = 1000)
        birch = cluster.Birch(n_clusters=3, threshold=1.7)
        sel_data=df_data[random_genes]

        kmeans.fit_predict(sel_data)
        birch.fit_predict(sel_data)

        # metrics - km & birch
        km_silhoute = metrics.silhouette_score(sel_data, kmeans.labels_, metric='euclidean')
        km_calinski_harabasz = metrics.calinski_harabasz_score(sel_data, kmeans.labels_)
        km_davies_bouldin = metrics.davies_bouldin_score(sel_data, kmeans.labels_)
        birch_silhoute = metrics.silhouette_score(sel_data, birch.labels_, metric='euclidean')
        birch_calinski_harabasz = metrics.calinski_harabasz_score(sel_data, birch.labels_)
        birch_davies_bouldin = metrics.davies_bouldin_score(sel_data, birch.labels_)

        # left for backwards compatibility
        if best_silhouette < birch_silhoute:
            best_genes = sel_data
            best_run = k
            best_silhouette = birch_silhoute

        # elbow_list.append(elbow_points)
        saved_config.append(sel_data)
        scores.append([kmeans.inertia_, km_silhoute, km_calinski_harabasz, km_davies_bouldin,
                        birch_silhoute, birch_calinski_harabasz, birch_davies_bouldin])
        stats.append([ sel_data.var(axis=1).values, sel_data.mean(axis=1).values,
                      sel_data.std(axis=1).values, sel_data.min(axis=1),
                      sel_data.max(axis=1).values])

    print("It took {:.5f} to run the experiment".format((time.time()-start_time_2)/60))

    # Create the df that we want to save
    scores_df = pd.DataFrame(scores, columns = ["inertia", "km_silhoute", "km_calinski","km_davies",
                                                "birch_silhoute", "birch_calinski","birch_davies"])
    stats_df = pd.DataFrame(stats, columns = ["var", "mean", "std", "min", "max"])

    results_all = scores_df.copy(deep=True)
    results_all["saved_config_genes"] = saved_config
    results_all = results_all.sort_values(by="birch_silhoute", ascending=False)
    results_dict = {
        "scores": scores_df,
        "all_scores": results_all,
        "best_genes": best_genes,
        "best_50_genes": results_all.iloc[:1000, :], #left the name for backwards compatibility
        "run_idx": best_run,
        "stats": stats_df }

    return results_dict

def process_random_exps(data, samples, config, plot_metrics=False, plot_tsne=False, individual_plots=False):
    """ Function that loads the random experiments objects and process them and by this we mean that the metrics plots, t-sne are created.

    Args:
        data (DataFrame): The dataset that's been used for experiments. We need this to computer the reference results, which in this case are represented by the most abundant genes
        samples (list): List of strings of the samples
        config (touple): Used for passing some config arguments like the colour pool or the base_path
        plot_metrics (bool, optional): If to plot the cluster metrics. Defaults to False.
        plot_tsne (bool, optional): If to plot the t-sne. Defaults to False.
        individual_plots (bool, optional): If True the t-SNE plots will be plot separately and not in subplots. Defaults to False.
    """

    df_metadata, df = data[0], data[1]
    base_path, colours_pool = config[0], config[1]
    df["mean"] = df.mean(axis=1)
    df["variance"] = df.var(axis=1)

    loaded_obj_1 = load_obj(base_path + "random_exp_" + str(100000))
    loaded_obj_2 = load_obj(base_path + "random_exp_" + str(100001))
    loaded_obj_3 = load_obj(base_path + "random_exp_" + str(100002))
    loaded_obj_4 = load_obj(base_path + "random_exp_" + str(100004))
    loaded_obj_5 = load_obj(base_path + "random_exp_" + str(100003))

    input_data = df.sort_values(by="variance", ascending=False).reset_index(drop=True).iloc[:197, :]
    input_data = input_data.drop(["mean", "variance"], axis=1)
    sel_results = run_random_config(input_data, 1, 197, samples)

    scores, best_genes, stats = [loaded_obj_1["scores"]], [loaded_obj_1["best_genes"]], [loaded_obj_1["stats"]]
    scores.append(loaded_obj_2["scores"])
    scores.append(loaded_obj_3["scores"])
    scores.append(loaded_obj_4["scores"])
    scores.append(loaded_obj_5["scores"])
    best_genes.append(loaded_obj_2["best_genes"])
    best_genes.append(loaded_obj_3["best_genes"])
    best_genes.append(loaded_obj_4["best_genes"])
    best_genes.append(loaded_obj_5["best_genes"])
    stats.append(loaded_obj_2["stats"])
    stats.append(loaded_obj_3["stats"])
    stats.append(loaded_obj_4["stats"])
    stats.append(loaded_obj_5["stats"])

    # find the common genes
    genes_indexes = [best_genes[0].columns.values]
    genes_indexes.append(best_genes[1].columns.values)
    genes_indexes.append(best_genes[2].columns.values)
    genes_indexes.append(best_genes[3].columns.values)
    genes_indexes.append(best_genes[4].columns.values)

    ref_genes = set(df.sort_values(by="variance", ascending=False).reset_index(drop=True).iloc[:197, :]["genes"])

    bst_gene_names = [set(df.iloc[genes_indexes[0]]["genes"])]
    bst_gene_names.append(set(df.iloc[genes_indexes[1]]["genes"]))
    bst_gene_names.append(set(df.iloc[genes_indexes[2]]["genes"]))
    bst_gene_names.append(set(df.iloc[genes_indexes[3]]["genes"]))
    bst_gene_names.append(set(df.iloc[genes_indexes[4]]["genes"]))

    common_genes = ref_genes & bst_gene_names[0] & bst_gene_names[1] & bst_gene_names[2] & bst_gene_names[3]

    print("Common genes to all: ", common_genes)
    print("Common genes to random: ", bst_gene_names[0] & bst_gene_names[1] & bst_gene_names[2])
    print("best random 1 - random 2: ", len(bst_gene_names[0] - bst_gene_names[1]))
    print("best random 2 - random 1: ", len(bst_gene_names[1] - bst_gene_names[0]))
    print("best random 1 - random 3: ", len(bst_gene_names[0] - bst_gene_names[2]))

    # # compare results
    best_genes_idxs = [loaded_obj_1["run_idx"]]
    best_genes_idxs.append(loaded_obj_2["run_idx"])
    best_genes_idxs.append(loaded_obj_3["run_idx"])
    best_genes_idxs.append(loaded_obj_4["run_idx"])
    best_genes_idxs.append(loaded_obj_5["run_idx"])

    # plot the overview of the metrig graphs. These contains the best configurations of each run
    dummy = pd.DataFrame( [scores[0].iloc[best_genes_idxs[0]].values,
                           scores[1].iloc[best_genes_idxs[1]].values,
                           scores[2].iloc[best_genes_idxs[2]].values,
                           scores[3].iloc[best_genes_idxs[3]].values,
                           scores[4].iloc[best_genes_idxs[4]].values,
                           sel_results["scores"].values[0]
                           ],
                        columns = scores[0].iloc[best_genes_idxs[0]].index.values)
    dummy[["km_silhoute", "birch_silhoute"]].plot(kind="bar", figsize=(17, 6))
    dummy[["km_calinski", "birch_calinski"]].plot(kind="bar", figsize=(17, 6))
    dummy[["km_davies", "birch_davies"]].plot(kind="bar", figsize=(17, 6))
    dummy["inertia"].plot(kind="bar", figsize=(17, 6))

    sel_inertia = sel_results["scores"].values[0][0]
    sel_silhouete_km = sel_results["scores"].values[0][1]
    sel_calinski_km = sel_results["scores"].values[0][2]
    sel_davies_km = sel_results["scores"].values[0][3]
    sel_silhouete_birch = sel_results["scores"].values[0][4]
    sel_calinski_birch = sel_results["scores"].values[0][5]
    sel_davies_birch = sel_results["scores"].values[0][6]

    # Plot the metrics scores
    if plot_metrics:
        for idx, score in enumerate(scores):
            score.reset_index(inplace=True)
            # inerita
            ax = score.plot(kind="scatter", y="inertia", x="index", figsize=(17, 6))
            ax.set_title("Inertia of Kmeans for Random experiment " + str(idx + 1))
            ax.plot(sel_inertia, 'ro')
            ax
            ax = score[["km_silhoute", "birch_silhoute"]].plot(kind="line", figsize=(17, 6))
            ax.set_title("KM and Birch, Silhoute scores for Random experiment " + str(idx + 1))
            ax.plot(sel_silhouete_birch, 'ro')
            ax.plot(sel_silhouete_km, 'bo')
            ax
            ax = score[["km_calinski", "birch_calinski"]].plot(kind="line", figsize=(17, 6))
            ax.set_title("KM and Birch, CalinsCalinski-Harabasz ki scores for Random experiment " + str(idx + 1))
            ax.plot(sel_calinski_birch, 'ro')
            ax.plot(sel_calinski_km, 'bo')
            ax
            ax = score[["km_davies", "birch_davies"]].plot(kind="line", figsize=(17, 6))
            ax.set_title("KM and Birch, Davies-Bouldin scores for Random experiment " + str(idx + 1))
            ax.plot(sel_davies_birch, 'ro')
            ax.plot(sel_davies_km, 'bo')
            ax
            # Silhouete
            ax = score.plot(kind="scatter", y=["km_silhoute"], x="index", figsize=(17, 6))
            ax.set_title("Kmeans Silhoute scores for Random experiment " + str(idx + 1))
            ax.plot(sel_silhouete_km, 'ro')
            ax
            ax = score.plot(kind="scatter", y=["birch_silhoute"], x="index", figsize=(17, 6))
            ax.set_title("birch, Silhoute scores for Random experiment " + str(idx + 1))
            ax.plot(sel_silhouete_birch, 'ro')
            ax
            # Calinski
            ax = score.plot(kind="scatter", y=["km_calinski"], x="index", figsize=(17, 6))
            ax.set_title("Kmeans Calinski scores for Random experiment " + str(idx + 1))
            ax.plot(sel_calinski_km, 'ro')
            ax
            ax = score.plot(kind="scatter", y=["birch_calinski"], x="index", figsize=(17, 6))
            ax.set_title("Birch Calinski scores for Random experiment " + str(idx + 1))
            ax.plot(sel_calinski_birch, 'ro')
            ax
            #Davies Bouie
            ax = score.plot(kind="scatter", y=["km_davies"], x="index", figsize=(17, 6))
            ax.set_title("Kmeans Davies scores for Random experiment " + str(idx + 1))
            ax.plot(sel_davies_km, 'ro')
            ax
            ax = score.plot(kind="scatter", y=["birch_davies"], x="index", figsize=(17, 6))
            ax.set_title("Birch Davies scores for Random experiment " + str(idx + 1))
            ax.plot(sel_davies_birch, 'ro')
            ax

    hover_data = ["samples", "diagnosis", "ter_avg", "ter_sd", "labels_cystometric_capacity","labels_date", "labels_ter_avg",  "labels_ter_sd"]
    if plot_tsne:
        best_genes_transp = []
        # Plot the PCA/t-sne
        birch_figs, kmeans_figs = [], []
        saved_data = pd.DataFrame(data=samples, columns=["samples"])
        for idx, bst_gene_name in enumerate(bst_gene_names):
            # finding the indexes equivalent in the original diff bu
            idxs = []
            for index, row in df.iterrows():
                if row["genes"] in list(bst_gene_name):
                    idxs.append(index)

            selected_genes = df.loc[idxs]
            selected_genes["genes"]= selected_genes["genes"]
            selected_genes["mean"] = selected_genes.mean(axis=1)
            selected_genes["variance"] = selected_genes.var(axis=1)
            best_genes_transp.append(selected_genes)
            config, percent = (samples, df_metadata, colours_pool), 1.0
            output_df, _ = selected_genes_experiment(selected_genes, percent, config, verbose=1, pca=False)

            # create the plot objects
            fig1 = px.scatter(output_df, x = "tsne_2d_one", y= "tsne_2d_two", color = "RawKMeans",
                              text = "samples", size="ter_avg",  hover_data = hover_data,
                              title = "K-means with t-SNE for Experiment run " + str(idx) )
            fig2 = px.scatter(output_df, x = "tsne_2d_one", y= "tsne_2d_two", color = "Birch", text = "samples",
                              size="ter_avg",  title = "Birch t-SNE for Experiment run " + str(idx),  hover_data = hover_data)
            fig1.update_traces(textposition='top right')
            fig2.update_traces(textposition='top right')

            kmeans_figs.append(fig1)
            birch_figs.append(fig2)
            saved_data["Birch_{}".format(idx)] = output_df["Birch"].values
            saved_data["Kmeans_{}".format(idx)] = output_df["Birch"].values

    _ = add_prcsd_metdata(saved_data, df_metadata)
    clusters_to_plot = saved_data.columns.values[1:len(bst_gene_names)*2+1]
    saved_data["tsne_2d_one"] = output_df["tsne_2d_one"]
    saved_data["tsne_2d_two"] = output_df["tsne_2d_two"]
    plot_clusters(saved_data, clusters_to_plot, "tsne_2d_one", "tsne_2d_two", colours_pool, "t-sne for random genes experiments", hover_data = hover_data) #marker_text="diagnosis"

    if individual_plots:
        for fig in birch_figs:
            fig.update_traces(textposition='top right')
            fig.show()
        for fig in kmeans_figs:
            fig.update_traces(textposition='top right')
            fig.show()

    # adding the genes names
    best_genes_transp = []
    for best_gene in best_genes:
        dummy = best_gene.transpose()
        dummy.columns = samples

# Process one single random exp
def process_single_exp(name, data, config, samples, show_figs = False, save_figs = False):
    """ Similar to process_random_exps, but it's for a single random experiment.

    Example:
        utilities.process_single_exp("random_exp_100001", data, config, diff_samples, show_figs= True, save_figs=False)

    Args:
        name (String): Name of  the file
        data (DataFrame): The dataset that's been used for experiments. We need this to computer the reference results, which in this case are represented by the most abundant genes
        config (touple): Used for passing some config arguments like the colour pool or the base_path
        samples (list): List of strings of the samples
        show_figs (bool, optional): If to display the figures. Defaults to False.
        save_figs (bool, optional): If to save the figures. Defaults to False.

    """
    df_metadata, df = data[0], data[1]
    base_path, colours_pool = config[0], config[1]
    loaded_obj = load_obj(base_path + name)

    df["mean"], df["variance"] = df.mean(axis=1), df.var(axis=1)

    input_data = df.sort_values(by="variance", ascending=False).reset_index(drop=True).iloc[:197, :]
    input_data = input_data.drop(["mean", "variance"], axis=1)
    sel_results = run_random_config(input_data, 1, 197, samples)

    scores, best_genes, best_genes_idxs = loaded_obj["scores"], loaded_obj["best_genes"], loaded_obj["run_idx"]
    genes_indexes = best_genes.columns.values
    ref_genes = set(df.sort_values(by="variance", ascending=False).reset_index(drop=True).iloc[:197, :]["genes"])

    # Check the common genes between reference and given experiment
    best_gene_names = set(df.iloc[genes_indexes]["genes"])
    print("Common genes to reference configurations: ", ref_genes & best_gene_names)

    dummy = pd.DataFrame( [scores.iloc[best_genes_idxs].values,
                           sel_results["scores"].values[0] ],
                        columns = scores.iloc[best_genes_idxs].index.values)
    dummy[["km_silhoute", "birch_silhoute"]].plot(kind="bar", figsize=(17, 6))
    dummy[["km_calinski", "birch_calinski"]].plot(kind="bar", figsize=(17, 6))
    dummy[["km_davies", "birch_davies"]].plot(kind="bar", figsize=(17, 6))
    dummy["inertia"].plot(kind="bar", figsize=(17, 6))

    #### Plot the PCA/t-sne
    # finding the indexes equivalent in the original diff bu
    idxs = []
    for index, row in df.iterrows():
        if row["genes"] in list(best_gene_names):
            idxs.append(index)
    selected_genes = df.loc[idxs]

    # apply PCA
    selected_genes["genes"]= selected_genes["genes"]
    selected_genes["mean"] = selected_genes.mean(axis=1)
    selected_genes["variance"] = selected_genes.var(axis=1)
    gen_exp_config, percent = (samples, df_metadata, colours_pool), 1.0
    output_df, _ = selected_genes_experiment(selected_genes, percent, gen_exp_config, verbose=1, pca=False)
    # t-sne plots
    hover_data = ["samples", "labels_cystometric_capacity","labels_date", "labels_ter_avg", "labels_ter_sd", "ter_avg", "ter_sd"]
    fig1 = px.scatter(output_df, x = "tsne_2d_one", y= "tsne_2d_two", color = "RawKMeans",
                      text = "samples", size="ter_avg", hover_data =  hover_data,
                      title = "K-means with t-SNE for Random Experiment" )
    fig2 = px.scatter(output_df, x = "tsne_2d_one", y= "tsne_2d_two", color = "Birch",
                      text = "samples", size="ter_avg", hover_data =  hover_data,
                      title = "Birch t-SNE for Random Experiment" )
    fig1.update_traces(textposition='top right')
    fig2.update_traces(textposition='top right')
    if show_figs:
        fig1.show()
        fig2.show()

    if save_figs:
        save_fig("random_exp/pca_kmeans_tsne_" + name, fig1)
        save_fig("random_exp/pca_birch_tsne_" + name, fig2)

# process the top performing genes, these are already sorted by birch_silhouete
def process_best_genes_config(file_path, df_bu, sel_samples):
    """ Process the data for the best gene configuration from the random experiments. This was created to quickly follow the steps. Note: this hasn't been update for more than 1 best configurations to be saved.

    Example of use: 
        utilities.process_best_genes_config(base_path+ "random_exp_" + str(100004), diff_bu_prcsd, diff_samples)

    Args:
        file_path (String): The file path of the experiments .csv
        df_bu ([DataFrame]): We need this to get the genes name. Make note this dataframe has to be unprocssed and unsorted with the except of eliminating some of the genes that have low expressed genes
        sel_samples ([String]): List of the samples.

    Returns:
        (DataFrame, list): A touple which contains the the scores in a dataframe and the second element is represented by the dataframe itself
    """
    results = load_obj(file_path)

    # this was kept for backwards compatibility, probably there are more than 50 configurations
    genes_configs = results["best_50_genes"]
    df_data = df_bu[sel_samples]
    sel_datas = []
    elbow_list = []
    for _, row in genes_configs.iterrows():
        # adding back the genes column
        gene_idxs = row["saved_config_genes"].columns.values
        sel_data = df_data.loc[gene_idxs]

        # calculate the elbow scores
        elbow_metrics = ["distortion", "silhouette", "calinski_harabasz"]
        elbow_points = []
        km = cluster.KMeans()
        for elbow_metric in elbow_metrics:
            visualizer = KElbowVisualizer(km, k=(2,12), metric=elbow_metric)
            visualizer.fit(sel_data)
            elbow_points.append(visualizer.elbow_value_)
            elbow_points.append(visualizer.elbow_score_)
        elbow_list.append(elbow_points)
        sel_data["genes"] = df_bu.loc[gene_idxs]["genes"]
        sel_datas.append(sel_data)

    elbow_df = pd.DataFrame(elbow_list, columns=["distortion_point", "distortion_score",
                                                 "silhouette_point", "silhouette_score",
                                                 "calinski_harabasz_point",
                                                 "calinski_harabasz_score"])

    return elbow_df, sel_datas

################ Saving objects ################
def save_fig(name, fig, width=None, height=None, scale=None, base_path=None):
    """ Saves a figure to the base path or othwerise to th specified path

    Args:
        name ([type]): Name of the figure
        fig ([type]): Figure object
        width ([type], optional): Default 1280.
        height ([type], optional): Default 720.
        scale ([type], optional): Default 2.
        base_path ([type], optional): Default root folder.
    """

    default_width, default_height, default_scale = 1280, 720, 2
    default_base_path =  "./data/figures"

    #override arguments if needed
    if width != None:
        default_width = width
    if height != None:
        default_height = height
    if scale != None:
        default_scale = scale
    if base_path != None:
        default_base_path = base_path

    fig.write_image(default_base_path + name + ".png", width=default_width, height=default_height, scale=default_scale)

def save_obj(filePath, obj):
    print("######## Saving object ########")
    pickle.dump(obj, open(filePath, "wb") )

def load_obj(filePath):
    if path.exists(filePath):
        loaded_obj = pickle.load(open(filePath, "rb"))
        print("######## Finished loading obj ########")
        return loaded_obj
    else:
        print("There was an error in loading the file")

def save_tissue_csv(df_bu_raw, df_meta, base_path, tissue = "D"):
    """ Function that saves the data for a specific tissue to a CSV file

    Args:
        df_bu_raw ([DataFrame]): The raw data (with all TPMs)
        df_meta ([DataFrame]): The metadata 
        base_path ([String]): Where to save it
        tissue (str, optional): [description]. Defaults to "D".
    """
    # get rid of the outlayer
    df_bu_raw.drop("Y2319-D", axis=1, inplace=True)

    # Filter the sample by the tissue type
    filtered_samples = []
    for col_name in df_bu_raw.columns[1:]:
        tissue_type = col_name.split("-")[1]
        if tissue == tissue_type:
            filtered_samples.append(col_name)

    # create the df with the filtered samples
    df_bu_filtered = pd.DataFrame()
    df_bu_filtered["samples"] = filtered_samples
    add_prcsd_metdata(df_bu_filtered, df_meta)

    filtered_samples.insert(0, "genes")
    df_bu_raw.loc[:, filtered_samples].to_csv(base_path + "_data.tsv", sep='\t', index=False)

    labelencoder= LabelEncoder()
    df_bu_filtered["batch"] = labelencoder.fit_transform(df_bu_filtered["diff_sequencer"])
    df_bu_filtered["dataset"] = df_bu_filtered["samples"]
    df_bu_filtered[["dataset", "batch"]].to_csv(base_path + "_bu_sequencer.tsv", sep='\t', index=False)

################ Others ################
def order_labels_size(df_series):
    """ This function ensures that the the cluster labels are always in the order of the cluster sizes. This means that the cluster label 0 will correspond for the largest cluster and for n-clusters the n-1 will be the cluster with the lowest members.

    Args:
        df_series (Pandas Series): The cluster labels to be mapped

    Returns:
        [Pandas Series]: New mappped pandas series
    """
    # ordered by the number of frequency
    cluster_count = [clust for clust, count in Counter(df_series).most_common()]
    new_labels = list(range(len(df_series.unique())))
    dic = dict(zip(cluster_count, new_labels))
    return df_series.map(dic)

def get_cluster_size(df_series):
    return [count for clust, count in Counter(df_series).most_common()]

def create_urls(genes):
    """ For a given list of genes strings generate the gene cards urls

    Args:
        genes ([String]): List of genes strings

    Returns:
        [String]: URls of gene cards
    """
    base_url = "https://www.genecards.org/cgi-bin/carddisp.pl?gene="
    urls =[]
    for gene in genes:
        urls.append(base_url + gene)

    return " ".join(urls)

def select_data(df, percent, new_idx_col = "genes"):
    no_genes = int(df.shape[0] * percent)
    input_data = drop_columns_if_needed(df.iloc[:no_genes])
    return input_data.set_index(new_idx_col)

# function taken from https://scikit-learn.org/stable/auto_examples/cluster/plot_agglomerative_dendrogram.html
def plot_dendrogram(model, **kwargs):
    # Create linkage matrix and then plot the dendrogram

    # create the counts of samples under each node
    counts = np.zeros(model.children_.shape[0])
    n_samples = len(model.labels_)
    for i, merge in enumerate(model.children_):
        current_count = 0
        for child_idx in merge:
            if child_idx < n_samples:
                current_count += 1  # leaf node
            else:
                current_count += counts[child_idx - n_samples]
        counts[i] = current_count

    linkage_matrix = np.column_stack([model.children_, model.distances_,
                                      counts]).astype(float)

    # Plot the corresponding dendrogram
    dendrogram(linkage_matrix, **kwargs)
