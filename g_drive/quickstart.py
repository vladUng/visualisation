#!/Users/vlad/opt/anaconda3/envs/visualisation python

from __future__ import print_function
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

# If modifying these scopes, delete the file token.pickle.
# https://developers.google.com/drive/api/v3/about-auth
# SCOPES = ['https://www.googleapis.com/auth/drive.metadata.readonly']
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
            '/keys/credentials.json', SCOPES)
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

def get_files(service, folder_id):

  # get the existing files
  existing_files = [f for f in os.listdir("saved_data/") if ".tsv" in f]

  # query to get the shared folder with data
  query="'%s' in parents" % (folder_id)
  response = service.files().list(q= query, spaces='drive',
                                  fields='nextPageToken, files(id, name)').execute()
  results = response.get('files', [])

  file_paths = []
  for result in results:
    file_name = str(result.get("name"))
    file_path = "saved_data/" + file_name

    if file_name not in existing_files:
      print("Need to download %s " % file_name)
      request = service.files().get_media(fileId=result.get("id"))
      fh = io.BytesIO()
      downloader = MediaIoBaseDownload(fh, request)
      done = False
      while done is False:
          status, done = downloader.next_chunk()
          print("Download %s %d%%." % (file_path, int(status.progress() * 100)))

      # saved files to the folder
      with open(file_path, "wb") as outfile:
        # Copy the BytesIO stream to the output file
        outfile.write(fh.getbuffer())

    file_paths.append(file_path)
  return file_paths

def process_files(file_paths):
  raw_dfs = {}
  for file_path in file_paths:
    file_name = file_path.split("/")[-1].split()[0]
    raw_dfs[file_name] = pd.read_csv(file_path, sep='\t')
  
  return raw_dfs

def main():
  # print("\n %s \n" % os.getcwd())
  drive_service = get_gdrive_access()
  folder_id = get_folder_id(drive_service, "test_folder")
  if folder_id != None:
    # donwload the files
    print("folder id", folder_id)
    file_paths = get_files(drive_service, folder_id)
    process_files(file_paths)

  # return the files paths


if __name__ == '__main__':
    main()