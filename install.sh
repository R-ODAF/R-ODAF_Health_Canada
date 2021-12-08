#!/usr/bin/env bash
# Find directory where script is running from
SCRIPT_DIR="$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"
# Copy default configuration to config.yaml
cp ${SCRIPT_DIR}/config/config.default.yaml ${SCRIPT_DIR}/config/config.yaml
# Modify the config file to update project root directory
sed -i 's:prefix.*$:prefix\: '${SCRIPT_DIR}'/:' config/config.yaml