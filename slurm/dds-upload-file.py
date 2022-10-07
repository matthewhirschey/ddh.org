import sys
import DukeDS

project = sys.argv[1]
filepath = sys.argv[2]

print("Uploading {} to {}".format(filepath, project))
DukeDS.upload_file(project, filepath)
