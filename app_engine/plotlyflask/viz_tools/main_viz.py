import numpy as np
import pandas as pd 
import plotlyflask.viz_tools.utilities as utilities
import plotly.express as px

from google.cloud import storage
import google.oauth2.id_token

from googleapiclient.discovery import build
from google_auth_oauthlib.flow import InstalledAppFlow
from googleapiclient.http import MediaIoBaseDownload
from google.auth.transport.requests import Request

import main
import os
import pickle
import io

from google_auth_oauthlib.flow import InstalledAppFlow
from google.auth.transport.requests import Request
import pickle
from googleapiclient.discovery import build
import os

SCOPES = ['https://www.googleapis.com/auth/drive', 
          'https://www.googleapis.com/auth/drive.file',
          'https://www.googleapis.com/auth/drive.metadata.readonly']


def get_gdrive_access():
  """
  Gets access to google drive
  """
  creds = None
  # The file token.pickle stores the user's access and refresh tokens, and is
  # created automatically when the authorization flow completes for the first
  # time.
  if os.path.exists('token.pickle'):
      with open('token.pickle', 'rb') as token:
          creds = pickle.load(token)

  # If there are no (valid) credentials available, let the user log in.
  if not creds or not creds.valid:
    if creds and creds.expired and creds.refresh_token:
        creds.refresh(Request())
    else:
        flow = InstalledAppFlow.from_client_secrets_file(
            'keys/client_id_desktop.json', SCOPES)
        creds = flow.run_local_server(port=0)
    # Save the credentials for the next run
    with open('token.pickle', 'wb') as token:
        pickle.dump(creds, token)
 
  service = build('drive', 'v3', credentials=creds)
  return service

def get_folder_id(service, folder_name):
    # query to get the shared folder with data
    query="sharedWithMe and name contains '%s'" % (folder_name)
    response = service.files().list(q= query, spaces='drive',
                                        fields='nextPageToken, files(id, name)').execute()
    results = response.get('files', [])
    if results:
      folder_id = results[0].get("id")
      return folder_id
    else:
      return None

def test_read_csv_gdrive(drive_service, to_download = []):
    dfs = [] #dataframes returned
    folder_id = get_folder_id(drive_service, "test_folder")   
    query="'%s' in parents" % (folder_id)

    response = drive_service.files().list(q= query, spaces='drive',
                                    fields='nextPageToken, files(id, name)', orderBy = "name").execute()
    results = response.get('files', [])

    for result in results:
        file_name = str(result.get("name"))

        if file_name in to_download:
            print("Need to download %s " % file_name)
            request = drive_service.files().get_media(fileId=result.get("id"))
            fh = io.BytesIO()
            downloader = MediaIoBaseDownload(fh, request)
            done = False
            while done is False:
                status, done = downloader.next_chunk()
                print("Downloaded %s %d%%." % (file_name, int(status.progress() * 100)))

            fh.seek(0)
            strings_io_obj = fh.read().decode('UTF-8')
            # df_metadata = pd.read_csv(strings_io_obj, sep='\t')
            dfs.append(pd.read_csv(io.StringIO(strings_io_obj), encoding='utf8', sep="\t"))

    return dfs

def test_bu_plots():
    #loading the files, filepaths and filenames
    base_path = "/Users/vlad/Documents/Code/York/BU_clustering/src/data/"

    bu_metadata = "200319_sequenced_sample_metadata.tsv" # metadata
    # Differentiate tissue type after batch
    some_diff_bc_raw = "some_diff_bc.tsv" #without Y2319 which hasn't differentiated properly (see BU_cluster_notebook)
    po_bc_raw = "p0_bc_2.tsv" #only p0

    drive_service = get_gdrive_access()
    to_download = [bu_metadata, some_diff_bc_raw, po_bc_raw]

    to_download.sort()
    dfs = test_read_csv_gdrive(drive_service, to_download)
    df_metadata, some_diff_bu_raw, p0_bu_raw = dfs[0], dfs[1], dfs[2]


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
