#!/bin/bash
set -euxo pipefail

echo "--- DEBUG INFO ---"
echo "Current PATH: $PATH"
echo "CONDA_PREFIX: $CONDA_PREFIX"
echo "--- END DEBUG INFO ---"

# This script runs after the reports environment is installed
# And adds the crosstools package, which isn't available through conda

mkdir -p output/.rodaf_internal/install_logs

# remove log file if it already exists
echo "Attempting to remove install_from_github.post-deploy.log"
rm -f output/.rodaf_internal/install_logs/install_from_github.post-deploy.log # Changed filename to match your redirect target.

# Execute Rscript directly from the CONDA_PREFIX bin directory.
# Check if Rscript exists before trying to run it
if [ ! -f "$CONDA_PREFIX/bin/Rscript" ]; then
    echo "ERROR: Rscript not found at $CONDA_PREFIX/bin/Rscript"
    ls -l "$CONDA_PREFIX"/bin/R* # List R related executables
    exit 1 # Exit with an error
fi

"$CONDA_PREFIX"/bin/Rscript workflow/envs/install_from_github.R >> output/.rodaf_internal/install_logs/install_from_github.post-deploy.log 2>&1
echo "Rscript command finished."

