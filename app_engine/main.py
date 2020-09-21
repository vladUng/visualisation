import datetime
# from flask import Flask, render_template, request, url_for, redirect, session
import flask 
import pickle
import os.path
import io
from httplib2 import Http
from googleapiclient.discovery import build
from google_auth_oauthlib.flow import InstalledAppFlow
from googleapiclient.http import MediaIoBaseDownload
from google.auth.transport.requests import Request
from apiclient import errors
import pandas as pd

from google.cloud import storage, datastore
import google.oauth2.id_token
from google.auth.transport import requests

import datetime




app = flask.Flask(__name__)
datastore_client = datastore.Client()
firebase_request_adapter = requests.Request()
app.secret_key = 'test'


# storing visits 
def store_time(dt):
    entity = datastore.Entity(key=datastore_client.key('visit'))
    entity.update({
        'timestamp': dt
    })

    datastore_client.put(entity)

def fetch_times(limit):
    query = datastore_client.query(kind='visit')
    query.order = ['-timestamp']

    times = query.fetch(limit=limit)

    return times

# buckets functions
def list_buckets():
    """Lists all buckets."""

    storage_client = storage.Client()
    buckets = storage_client.list_buckets()

    for bucket in buckets:
        print(bucket.name)

def upload_blob(source_file_name, destination_blob_name):
    """Uploads a file to the bucket."""
    bucket_name = "visualisation-jbu.appspot.com"
    # source_file_name = "local/path/to/file"
    # destination_blob_name = "storage-object-name"

    storage_client = storage.Client()
    bucket = storage_client.bucket(bucket_name)
    blob = bucket.blob(destination_blob_name)

    blob.upload_from_file(source_file_name)

    print(
        "File {} uploaded to {}.".format(
            source_file_name, destination_blob_name
        )
    )

# Google drive functions
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
  # pickle_path = '/tmp/token.pickle' #for
  #  app engine
  pickle_path = 'token.pickle'
  if os.path.exists(pickle_path):
      with open(pickle_path, 'rb') as token:
          creds = pickle.load(token)

  # If there are no (valid) credentials available, let the user log in.
  if not creds or not creds.valid:
    if creds and creds.expired and creds.refresh_token:
        creds.refresh(Request())
    else:
        flow = InstalledAppFlow.from_client_secrets_file(
            'keys/client_id_web.json', SCOPES)

        # flow.redirect_uri = url_for("auth", _external=True)
        # flow.redirect_uri = flask.url_for("oauth2callback", _external=False)
        flow.redirect_uri = 'https://visualisation-jbu.nw.r.appspot.com/oauth2callback'
        # creds = flow.run_local_server(port=0)

        # Generate URL for request to Google's OAuth 2.0 server.
        # Use kwargs to set optional request parameters.
        authorization_url, state = flow.authorization_url(
            # Enable offline access so that you can refresh an access token without
            # re-prompting the user for permission. Recommended for web server apps.
            access_type='offline',
            # Enable incremental authorization. Recommended as a best practice.
            include_granted_scopes='true')

        # Store the state so the callback can verify the auth server response.
        flask.session['state'] = state
        print("authroisation url", authorization_url)
        print("State: ", state)
        return flask.redirect("www.google.com")
 
    # Save the credentials for the next run
    with open(pickle_path, 'wb') as token:
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

def get_files(service, folder_id):

  # get the existing files
  existing_files = get_existing_files()

  # query to get the shared folder with data
  query="'%s' in parents" % (folder_id)
  response = service.files().list(q= query, spaces='drive',
                                  fields='nextPageToken, files(id, name)').execute()
  results = response.get('files', [])

  new_files_added = []
  for result in results:
    file_name = str(result.get("name"))

    if file_name not in existing_files:
      print("Need to download %s " % file_name)
      request = service.files().get_media(fileId=result.get("id"))
      fh = io.BytesIO()
      downloader = MediaIoBaseDownload(fh, request)
      done = False
      while done is False:
          status, done = downloader.next_chunk()
          print("Downloaded %s %d%%." % (file_name, int(status.progress() * 100)))

      fh.seek(0)
      upload_blob(fh, file_name)
      new_files_added.append(file_name)

  return new_files_added

def get_existing_files():

  # for bucket list 
  storage_client = storage.Client()
  blobs = storage_client.list_blobs("visualisation-jbu.appspot.com", delimiter="/")

  filenames = []
  for blob in blobs:
     filenames.append(blob.name)

  return filenames
  # for localhost
  # return [f for f in os.listdir("../saved_data/") if ".tsv" in f]

def process_files(file_paths):
  raw_dfs = {}
  for file_path in file_paths:
    file_name = file_path.split("/")[-1].split()[0]
    raw_dfs[file_name] = pd.read_csv(file_path, sep='\t') 
  
  return raw_dfs


@app.route('/')
def root():
  # pickle_path = '/tmp/token.pickle' #for app engine
  pickle_path = 'token.pickle'
  if os.path.exists(pickle_path):
    with open(pickle_path, 'rb') as token:
        creds = pickle.load(token)
    print("Alread having the credentials")
    flask.session['credentials'] = {
      'token': creds.token,
      'refresh_token': creds.refresh_token,
      'token_uri': creds.token_uri,
      'client_id': creds.client_id,
      'client_secret': creds.client_secret,
      'scopes': creds.scopes}
    return flask.redirect(flask.url_for("test_api_request"))
  else:
    return flask.redirect(flask.url_for("oauth2"))


  return flask.render_template('templates/index.html', times=[1,2,3])


  # Google sign in implementation
  id_token = flask.request.cookies.get("token")
  error_message = None
  claims = None
  times = [None]
  if id_token:
    try:
        # Verify the token against the Firebase Auth API. This example
        # verifies the token on each page load. For improved performance,
        # some applications may wish to cache results in an encrypted
        # session store (see for instance
        # http://flask.pocoo.org/docs/1.0/quickstart/#sessions).
        claims = google.oauth2.id_token.verify_firebase_token(
            id_token, firebase_request_adapter)
    except ValueError as exc:
        # This will be raised if the token is expired or any other
        # verification checks fail.
        error_message = str(exc)

    # Record and fetch the recent times a logged-in user has accessed
    # the site. This is currently shared amongst all users, but will be
    # individualized in a following step.
    store_time(datetime.datetime.now())
    times = fetch_times(10)
  print(claims, error_message, times)
  return flask.render_template(
        'index_2.html',
        user_data=claims, error_message=error_message, times=times)


@app.route('/oauth2')
def oauth2():

  print("\n\n#Start the auth process")
  # Create flow instance to manage the OAuth 2.0 Authorization Grant Flow steps.
  flow = flow = InstalledAppFlow.from_client_secrets_file('keys/client_id_web.json', SCOPES)

  # The URI created here must exactly match one of the authorized redirect URIs
  # for the OAuth 2.0 client, which you configured in the API Console. If this
  # value doesn't match an authorized URI, you will get a 'redirect_uri_mismatch'
  # error.
  # flow.redirect_uri = flask.url_for('oauth2callback', _external=True)
  # flow.redirect_uri = 'https://visualisation-jbu.nw.r.appspot.com/oauth2callback'
  flow.redirect_uri = flask.url_for('oauth2callback', _external=True)

  authorization_url, state = flow.authorization_url(
      # Enable offline access so that you can refresh an access token without
      # re-prompting the user for permission. Recommended for web server apps.
      access_type='offline',
      # Enable incremental authorization. Recommended as a best practice.
      include_granted_scopes='true')

  # Store the state so the callback can verify the auth server response.
  flask.session['state'] = state
  return flask.redirect(authorization_url)

@app.route('/oauth2callback')
def oauth2callback():
  print("\n\n Oauth callback")
  print("State: ", flask.session)
  # if flask.session['state']:
  state = flask.session['state']

  flow = InstalledAppFlow.from_client_secrets_file(
            'keys/client_id_web.json', SCOPES, state=state)
  # flow = InstalledAppFlow.from_client_secrets_file(
  #     'client_secret.json',
  #     scopes=['https://www.googleapis.com/auth/drive.metadata.readonly'],
  #     state=state)
  # flow.redirect_uri = flask.url_for('oauth2callback', _external=True)

  flow.redirect_uri = flask.url_for('oauth2callback', _external=True)
  # Use the authorization server's response to fetch the OAuth 2.0 tokens.
  authorization_response = flask.request.url

  # authorization_response = "https://visualisation-jbu.nw.r.appspot.com/oauth2callback"
  flow.fetch_token(authorization_response=authorization_response)

  # Store the credentials in the session.
  # ACTION ITEM for developers:
  #     Store user's access and refresh tokens in your data store if
  #     incorporating this code into your real app.
  credentials = flow.credentials
  flask.session['credentials'] = {
      'token': credentials.token,
      'refresh_token': credentials.refresh_token,
      'token_uri': credentials.token_uri,
      'client_id': credentials.client_id,
      'client_secret': credentials.client_secret,
      'scopes': credentials.scopes}
  print("\n##Auth finished")

  # Save the credentials for the next run
  # pickle_path = '/tmp/token.pickle' #for app engine
  pickle_path = 'token.pickle'
  with open(pickle_path, 'wb') as token:
    pickle.dump(credentials, token)

  return flask.redirect(flask.url_for('test_api_request'))


@app.route('/test_api_request')
def test_api_request():
  print("\n\n#####Testing api call: ")
  if 'credentials' not in flask.session:
    return flask.redirect('authorize')
# 
  # Load credentials from the session.
  credentials = google.oauth2.credentials.Credentials(
      **flask.session['credentials'])

  drive_service = build('drive', 'v3', credentials=credentials)

  folder_id = get_folder_id(drive_service, "test_folder")
  if folder_id != None:
    # donwload the files
    print("folder id", folder_id)
    files_added = get_files(drive_service, folder_id)
    print("\n\n##Files added: ", files_added)
  
  # Save credentials back to session in case access token was refreshed.
  # ACTION ITEM: In a production app, you likely want to save these
  #              credentials in a persistent database instead.
  flask.session['credentials'] = credentials_to_dict(credentials)

  return flask.render_template('templates/index.html', times=files_added)


def credentials_to_dict(credentials):
  return {'token': credentials.token,
          'refresh_token': credentials.refresh_token,
          'token_uri': credentials.token_uri,
          'client_id': credentials.client_id,
          'client_secret': credentials.client_secret,
          'scopes': credentials.scopes}
          
if __name__ == '__main__':
    # This is used when running locally only. When deploying to Google App
    # Engine, a webserver process such as Gunicorn will serve the app. This
    # can be configured by adding an `entrypoint` to app.yaml.
    # Flask's development server will automatically serve static files in
    # the "static" directory. See:
    # http://flask.pocoo.org/docs/1.0/quickstart/#static-files. Once deployed,
    # App Engine itself will serve those files as configured in app.yaml.

    # When running locally, disable OAuthlib's HTTPs verification.
    # ACTION ITEM for developers:
    #     When running in production *do not* leave this option enabled.
    os.environ['OAUTHLIB_INSECURE_TRANSPORT'] = '1'

    app.run(host='127.0.0.1', port=8080, debug=True)