depmap-shiny
============

Openshift deployment of Depmap application using shiny-server

## Contents

- CreateGeneSummaryJob.yaml: Kubernetes batch job to create `gene_summary.RData` - WIP
- application.yaml: Tweaked output of `process-template.sh`
- process-template.sh: Script using [juhahu/shiny-openshift](https://github.com/juhahu/shiny-openshift/) to bootstrap Openshift resources in `application.yaml`
- params: Parameters used in the template, named in `process-template.sh`

## Deployment

1. Create Openshift objects (BuildConfig, DeploymentConfig, Route, Service, PersistentVolumeClaim): `oc create -f application.yaml`
2. Wait for build to complete and deployment to start
3. Place data files on the PersistentVolume: e.g. `oc cp ./data/gene_summary.RData depmap-shiny-1-xxxxx:/srv/data/`
4. Visit the public route at http://depmap-shiny.apps.cloud.duke.edu to verify application is running

## Notes / TODOs

- The template is a good way to get shiny server up and running, but includes examples and packages we do not need
- The shiny application is built from Dockerfile.shiny, but the template does not allow for customization of the `dockerfilePath`, so this has been manually added to the BuildConfig in `application.yaml`.
- The CreateGeneSummaryJob does start and uses the NCBI API key (it must be in a secret), but in my testing on dev, the job did not finish, so instead I've just provided the
