from ddsc.sdk.client import Client
from ddsc.sdk.client import PathToFiles
import json

client = Client()
project = client.get_project_by_name('ddh-data')
ptf = PathToFiles()
ptf.add_paths_for_children_of_node(project)

files = []

for p, f in ptf.paths.items():
  files.append({'dest': p, 'source': f.id, 'type': 'DukeDS'})

print(json.dumps({"items": files}))
