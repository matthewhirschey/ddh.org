import DukeDS
import os

quarter = os.getenv('DDH_QUARTER', '20Q2')

ddh_files = DukeDS.list_files('ddh-com-data')
for ddh_file in ddh_files:
  if ddh_file.startswith("/data/{}".format(quarter)):
    print(ddh_file)
    DukeDS.download_file('ddh-com-data', ddh_file)
