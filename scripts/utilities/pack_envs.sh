#!/bin/bash

# Use jq to parse the JSON and extract the environment directories
ENV_PATHS=$( conda env list --json | jq -r '.envs[] ' | grep R-ODAF)

# add a line to print out every member of ENV_PATHS
echo $ENV_PATHS >&2


# Loop through the paths and pack each environment
for ENV_PATH in $ENV_PATHS; do
    echo $ENV_PATH >&2
    # Use the basename of the environment path as the name of the packed archive
    ENV_NAME=$(basename $ENV_PATH)
    # Pack the environment
    conda-pack --ignore-missing-files -p $ENV_PATH -o ./env_pack.${ENV_NAME}.tar.gz
done