#!env bash
set -euxo pipefail

# This script runs after the reports environment is installed
# And adds the crosstools package, which isn't available through conda

mkdir -p output/.rodaf_internal/install_logs

# remove log file if it already exists
rm -f output/.rodaf_internal/install_logs/install_from_github.post-deploy.log # Changed filename to match your redirect target.

# Check if Rscript exists in CONDA_PREFIX binbefore trying to run it
if [ ! -f "$CONDA_PREFIX/bin/Rscript" ]; then
    echo "ERROR: Rscript not found at $CONDA_PREFIX/bin/Rscript"
    ls -l "$CONDA_PREFIX"/bin/R* # List R related executables
    exit 1 # Exit with an error
fi

# Run install script with logging
"$CONDA_PREFIX"/bin/Rscript workflow/envs/install_from_github.R >> output/.rodaf_internal/install_logs/install_from_github.post-deploy.log 2>&1

echo "Rscript command finished."

