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

# If modifying these scopes, delete the file token.pickle.
# https://developers.google.com/drive/api/v3/about-auth
# SCOPES = ['https://www.googleapis.com/auth/drive.metadata.readonly']
SCOPES = ['https://www.googleapis.com/auth/drive',
          'https://www.googleapis.com/auth/drive.file',
          'https://www.googleapis.com/auth/drive.metadata.readonly']

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


  # query to get the shared folder with data
  query="'%s' in parents" % (folder_id)
  response = service.files().list(q= query, spaces='drive',
                                  fields='nextPageToken, files(id, name)').execute()
  results = response.get('files', [])

  for result in results:
    file_id = result.get('id')
    file_name = result.get("name")
    print(file_name)

    request = service.files().get_media(fileId=file_id)
    fh = io.BytesIO()
    downloader = MediaIoBaseDownload(fh, request)
    done = False
    while done is False:
        status, done = downloader.next_chunk()
        print("Download %d%%." % int(status.progress() * 100))

    # saved files to the folder
    with open("/Users/vlad/Documents/Code/York/visualisation/visualisation/saved_data/" + file_name, "wb") as outfile:
      # Copy the BytesIO stream to the output file
      outfile.write(fh.getbuffer())


def main():
  """Shows basic usage of the Drive v3 API.
  Prints the names and ids of the first 10 files the user has access to.
  """
  creds = None
  # The file token.pickle stores the user's access and refresh tokens, and is
  # created automatically when the authorization flow completes for the first
  # time.
  if os.path.exists('visualisation/token.pickle'):
      with open('visualisation/token.pickle', 'rb') as token:
          creds = pickle.load(token)
  # If there are no (valid) credentials available, let the user log in.
  if not creds or not creds.valid:
      if creds and creds.expired and creds.refresh_token:
          creds.refresh(Request())
      else:
          flow = InstalledAppFlow.from_client_secrets_file(
              'visualisation/keys/credentials.json', SCOPES)
          creds = flow.run_local_server(port=0)
      # Save the credentials for the next run
      with open('visualisation/token.pickle', 'wb') as token:
          pickle.dump(creds, token)

  service = build('drive', 'v3', credentials=creds)

  folder_id = get_folder_id(service, "test_folder")
  if folder_id != None:
    # donwload file
    print("folder id", folder_id)
    get_files(service, folder_id)


if __name__ == '__main__':
    main()