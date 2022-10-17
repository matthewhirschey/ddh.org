#!/bin/bash
export AWS_DEFAULT_PROFILE=ddh
#copilot init --app ddh-org \
#  --name ddh-org \
#  --type "Load Balanced Web Service" \
#  --dockerfile "./Dockerfile.simple" \
#  --deploy
copilot deploy --app ddh-org
