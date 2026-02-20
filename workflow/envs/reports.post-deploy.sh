#!env bash
set -o pipefail

# This script runs after the reports environment is installed
# And adds the crosstools package, which isn't available through conda

mkdir -p output/.rodaf_internal/install_logs

# remove log file if it already exists
rm -f output/.rodaf_internal/install_logs/install_from_github.post-deploy.log # Changed filename to match your redirect target.

# Execute Rscript directly from the CONDA_PREFIX bin directory.
"$CONDA_PREFIX"/bin/Rscript workflow/envs/install_from_github.R >> output/.rodaf_internal/install_logs/install_from_github.post-deploy.log 2>&1