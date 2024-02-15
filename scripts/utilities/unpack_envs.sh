#!/bin/bash

# Directory where the packed environments will be unpacked
UNPACK_DIR="./.snakemake/conda"

# Create the directory if it doesn't exist
mkdir -p $UNPACK_DIR

# Loop over the tar.gz files and unpack each environment
for ENV_ARCHIVE in ./env_pack*.tar.gz; do
    # Extract the environment name from the archive name
    ENV_NAME=$(basename $ENV_ARCHIVE .tar.gz | sed 's/env_pack//')

    # Create a directory for the environment
    ENV_DIR="${UNPACK_DIR}/${ENV_NAME}"
    mkdir -p $ENV_DIR

    # Unpack the environment
    tar -xzf $ENV_ARCHIVE -C $ENV_DIR

    # Activate the unpack script provided by conda-pack
    source ${ENV_DIR}/bin/activate
    conda-unpack

    # Clean up the archive file
    rm $ENV_ARCHIVE
done