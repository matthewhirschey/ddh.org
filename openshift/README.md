ddh-app
=======

Openshift deployment of ddh application using shiny server

## Contents

- Build.yaml: Build configurations and image streams for base image and application image
- Deployment.yaml: Deployment configuration to run application pods
- DownloadDDHDataJob.yaml: Batch Job for downloading data from DDS
- Storage.yaml: Configuration for requesting persistent storage to house data
- file-list.json: List of files to fetch from DDS
- make-file-list.py: Script to produce `file-list.json` using DukeDSClient

## Deployment

1. Create shared storage `oc create -f Storage.yaml`
2. Download data sets `oc create -f DownloadData.yaml` and wait for job to complete
3. Create Build and image configurations `oc create -f Build.yaml`
4. Deploy application: `oc create -f Deployment.yaml`
5. Ask openshift for the application route: `oc get route ddh-shiny-route`

## Development

The build/deployment is organized into two Docker images:

1. A **base** image with R, shiny-server, and packages installed: [base-shiny/Dockerfile](../base-shiny/Dockerfile)
2. An **application** image that builds _from_ the base but just adds the ddh code: [Dockerfile.ddh-shiny-app](../Dockerfile.ddh-shiny-app)

The base image takes about 25 minutes to build, and the application image builds much more quickly (just a few seconds). So we've configured code changes to just rebuild the application image to speed up the deployment process.

If additional Linux or R packages are needed, they should be added to the base image, since they are unlikely to change very often.

Note that the application image build is linked to the base image, so rebuilding the base will trigger a rebuild of the application.
