#!env bash
set -o pipefail

# This script runs after the reports environment is installed
# And adds the crosstools package, which isn't available through conda
conda activate -p $CONDA_PREFIX

mkdir -p output/.rodaf_internal/install_logs

(test -f output/.rodaf_internal/install_logs/myRscript.post-deploy.log || rm output/.rodaf_internal/install_logs/myRscript.post-deploy.log)

$CONDA_PREFIX/bin/Rscript workflow/envs/install_from_github.R >> output/.rodaf_internal/install_logs/install_from_github.post-deploy.log 2>&1