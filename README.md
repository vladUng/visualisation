# 1. JBU Visualisation tool

This is a tool that it is used by the JBU unit to visualise a gene's TPM's values across different datasets. The application runs at: `http://localhost:8080/dashapp/`

<!-- TOC -->

- [1. JBU Visualisation tool](#1-jbu-visualisation-tool)
- [2. Requirements](#2-requirements)
- [3. Installation](#3-installation)
  - [3.1. Adding a new dataset](#31-adding-a-new-dataset)
- [4. Features](#4-features)
- [5. Wishlist](#5-wishlist)
- [6. Notes on implementation](#6-notes-on-implementation)
  - [6.1. Technologies](#61-technologies)
  - [6.2. Why Flask + Dash](#62-why-flask--dash)
  - [6.3. Architecture](#63-architecture)
  - [6.4. Tips](#64-tips)
- [7. (Google) Cloud integration](#7-google-cloud-integration)

<!-- /TOC -->

# 2. Requirements

* Python 3.7. Specifically for Windows 64bit version
* Python packages:
  * Flask 
  * python-dotenev
  * dash
  * plotly
  * psutil
  * pandas

# 3. Installation 

1. Download the .tsv data from whereve it is stored (at the writting of this, was stored on Google Drive)
2. Go to the root folder of the project and run the following `pip install -r requirements.tx`. This will install all the python dependencies 
3. After everything is installed:
   1. Linux/Mac: `python wsgi.py`
   1. Windows: double click `wsgi.py`

## 3.1. Adding a new dataset

At the moment of writing there are 11 .tsv files, stored in `/data` folder and can be split in a few categories:
* A metadata file,
* A master .tsv file, which acts like a record file to keep track of the datasets used. It's also used to create an all data file: "all_data.tsv"
*  A all_data.tsv file
*  Rest of the 8 .tsv files are various datasets that hold the TPM values for each of the sample.

When a new dataset is added:
1. Download the new .tsv file to `/data`
2. Update the master.tsv with the name of the new .tsv
3. Run the the Python function `update_all_datasets("data/"`. This will update all_data.tsv file with the new TPMs / Samples

# 4. Features

* Do an exact search on the gene entered and plot the TPMs value 
* When one dataset is select, it plots the shared number column values
* Change the x-axis select between: Sample, shared_num_same_col, Tissue, Dataset,NHU_differentiation, Gender, TER, Substrate
* Filter out the loose TER barrier 
* Export the plot to PDF and that data to .CSV
* Select different type of figure plots: Swarm, violin, violin and points, box, box and points
* Choose between linear and log10 axis 

# 5. Wishlist 

Below are features to be implemented: 

* mann whitney, dataset pairwise comparison

# 6. Notes on implementation

## 6.1. Technologies 
The application is build on [Dash](https://plotly.com/dash/) which is the client-server solution for the interactive plotting library [Plotly](https://plotly.com/). Esentially, Dash is a wrapper arround the [Flask](https://flask.palletsprojects.com/en/1.1.x/) micro-server framework (if you are familiar with Node.js, it's similar to Restify). This means that on top of the Plotly graphing tools we can now create a GUI on a web-browser. See the [documentation](https://dash.plotly.com/layout) on Dash about how to add elements (it's actually good).

## 6.2. Why Flask + Dash
Even though Dash is a wrapper around a Flask (i.e. that a Flask server is started with a Dash app) and there is not an explicit need to run the Flask, we do that. The reason for this is that we initially wanted to integrate the tool in Google Cloud Service (more on this Cloud section) and we needed a little bit more freedom than the Dash let us. 

## 6.3. Architecture 
I have followed [this tutorial](https://hackersandslackers.com/plotly-dash-with-flask/) on how to do the integration. You don't need to follow the blog if you don't want to get the details, the important points are below:
* The `wsgi.py` does the initialisation and starts the server. This includes running the init function from the root folder. 
* The init script runs the flask server with the configuration from `config.py` that tells flask where to look for the resources (e.g. CSS files, html, images etc.). Also, it starts the Dash app
* The `/plotlyflask` folder contains all the code to make the visualisation tool to work.
* From an UI persepctive the tool is divided in the following: Gene search field, Dataset panel, Metadata panel, figure and export panel. The panels are constructing in each of their function and can be modified independently of the other elements.
* All the interactiactions with the html elements (i.e. button pressing, ticking etc.) are handled in init_callbacks function. Note, this is a workaround of the normal Dash implementation with Flask, normally you wouldn't have the callbacks inside a function, but outside and the dash app object global. 
* When there is nothing to plot, for situations when the gene is not found or at the start of the tool we return to the figure output `{ "data" : {} }`. The reason for this is that the plotly figures are dictionaries and we just return a figure with no data. This helps to deal with the cases when we want to manipulate the figure, like exporting. 
* **Importantly!** As we want the columns with _incl values as datasets (this means that we need to find the samples that have that column value = "Y"), we simply duplicate the samples and added to the original data. This involved to duplicate both the samples in metadata and in the `all_data.tsv`. We've marked the samples that are duplicated by the following rule "_incl_" + "including column initials". 
  * When we export the data to .csv we removed the duplicated samples 
* The two plotting functions are very similar and probably can be abastrated further, but we wanted to keep seperately when we extend the functionality for the metadata panel
* The `viz_tools/` folder was intended to integrate some of the functionality from the BU clustering project to this one, but we haven't actually completed that. Moreover, not sure if it's required, but left the files just in case that will be needed. 

## 6.4. Tips

* If you create a new file and need access to the flask app, simply import `from flask import current_app as app`
* An html element value can be se only by one callback (i.e. Output on the dash_app decorators) 


# 7. (Google) Cloud integration

Below are some of my notes from when I've tried to integrate the visualisation tool with Google Cloud. Some of them may be a bit outadate so it's better to consult the Google's oficial documentation. Also, the git repo's history may be useful. We've used the following services (which you need to set up first on your Google Cloud Project -GCP-):
* Google App Engine
* Google Storage Cloud
* Google Drive API
* Google Sign-in

At the begining of the project we've tried to create the tool directly on the Google Cloud. We aimed to that with the following constraints: efficient, fast, very cheap (close to free), secure. Unfortunatelly, we couldn't find a solution that satisfed all of these in the given time framed. 

We have tried to host the server on serverless solution called Google App Engine which basically runs a script when a request to its instance address is called. The server won't be always on and will restart automatically if needed. Some important aspects about Google App Engine
* The default resources alocated are small, for example you have a few MBs for RAM and the CPU is in the range of MHz. This worked fine for our current data but it may represented an issue in the feature 
* The data can't be stored and if so it will be temporary on the RAM which is not very large.
  * To solve this we have integrated Google Storage Cloud which is used to store data. As this is another google service has the potential to add extra costs. 
  * The initial plan was to host the data on Google Drive and move it to Google Storage Cloud. We have chosen this, as there is an API for both Google Services and the Google Drive is known by all the university members.
  * We've also tried to stream directly a .tsv files directly to a DataFrame which is possible as the `pd.read_csv` receives either a filepath or a StringIO. The google drive API returns a BufferIO which can be converted to StringIO.  <-- This may be relevant for future 
* Google Sign-in, Google drive authorisation integration - we needed to do this in order to insure security.

Google Drive API it's very easy to use but the documentation is poor in `g_drive/quickstart.py` is the integration of some basic api calls like:
* authorise the google drive
* getting a folder_id of a shared folder from Google Drive
* get the files from a folder_id
* save those files

In the `g_drive/main.py` there are some additional features such as:
* integration with Google Storage Cloud for:
  * uploading data
  * get the files from the buckets