#!/usr/bin/env bash

oc process \
  -f https://raw.githubusercontent.com/juhahu/shiny-openshift/master/shinytemplate.yml \
  --param-file params \
  -o yaml

