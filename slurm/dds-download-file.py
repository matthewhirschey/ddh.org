import sys
import DukeDS

project = sys.argv[1]
remotepath = sys.argv[2]
filepath = sys.argv[3]

print("Downloading {}/{} to {}".format(project, remotepath, filepath))
DukeDS.download_file(project, remotepath, filepath)
