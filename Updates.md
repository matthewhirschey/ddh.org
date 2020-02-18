# Steps to apply quarterly update

## Initial code changes
- Update release name, urls, and na_cutoff in [code/current_release.R]( https://github.com/hirscheylab/ddh/blob/master/code/current_release.R)
- Update [code/methods.md](https://github.com/hirscheylab/ddh/blob/master/code/methods.md) and supporting images
- Update `DMVER` in Makefile [code/Makefile](https://github.com/hirscheylab/ddh/blob/98e53f956439c570aefa7b8b2583cee1f84e8b2e/Makefile#L13) to match the new release name
- PR and Merge the above changes

NOTE: If the dependencies used for data generation have changed you will need to wait for the image to [finish buliding on dockerhub](https://hub.docker.com/repository/docker/dukegcb/ddh).

## Generate Data
All these steps are done on a Slurm cluster.
You will need at least 4G of storage run these steps.

Clone this repo, change into the directory and create a directory to hold slurm logs:
```
   git clone git@github.com:hirscheylab/ddh.git
   cd ddh   
   mkdir logs
```
Create a file named `config.sh` with the following contents (replacing TODO_NCBI_KEY with a valid key):
```
export DDH_BASE="."
SINGULARITY_BASE="${DDH_BASE}/singularity"
export SINGULARITY_TMPDIR="${SINGULARITY_BASE}/tmp"
export SINGULARITY_CACHEDIR="${SINGULARITY_BASE}/cache"
export SINGULARITY_IMAGEDIR="${SINGULARITY_BASE}/images"
export ENTREZ_KEY="TODO_NCBI_KEY"
```

Run build slurm job:
```
sbatch build-slurm.sh
```

After about 10 hours it should have finished.
View files in `logs` directory.
Re-run build-slurm.sh if it failed.

After successful upload data to DukeDS:
```
sbatch upload-slurm.sh
```

You can now delete the ddh directory from the slurm cluster.

## Deploy Data
To deploy the files requires access to the openshift okd console and a terminal environment with the [openshift oc command installed](https://docs.okd.io/latest/cli_reference/get_started_cli.html).

They should be done from within the `openshift` directory of a clone of this repository.

Update the [list of files to upload](https://github.com/hirscheylab/ddh/blob/master/openshift/file-list.json) by running:
```
python make-file-list.py > file-list.json
```

Use the `Copy login command` option in the top right corner of the openshift okd console.
Open a terminal and paste the copied command to setup `oc` to connect to the openshift cluster.

Switch `oc` to point at the depmap project.
```
oc project depmap
```

Delete the old list of files from openshift depmap.
```
oc delete configmap ddh-data-file-list
```

Create the new list of files from openshift depmap.
```
oc create configmap ddh-data-file-list --from-file=file-list.json
```
TODO: Test the above command to see if it works

Delete the old stage data job:
```
oc delete job download-ddh-data
```

Create a new stage data job:
```
oc create -f DownloadJob.yaml 
```

## Save file-list.json
Commit/PR/Merge changes to `file-list.json`.



