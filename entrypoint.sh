#!/usr/bin/env bash
set -e

if [[ -z "$HOST_UID" ]]; then
    echo "ERROR: please set HOST_UID" >&2
    exit 1
fi
if [[ -z "$HOST_GID" ]]; then
    echo "ERROR: please set HOST_GID" >&2
    exit 1
fi

# Use this code if you want to modify an existing user account:
groupmod --gid "$HOST_GID" R-ODAF
usermod --uid "$HOST_UID" R-ODAF
# Necessary to ensure all files are updated to the new UID and GID
chown -R R-ODAF:R-ODAF /home/R-ODAF/R-ODAF_Health_Canada

# Drop privileges and execute next container command, or 'bash' if not specified.
if [[ $# -gt 0 ]]; then
    exec gosu R-ODAF "$@"
    #exec sudo -u -H R-ODAF -- "$@"
else
    exec gosu R-ODAF "bash"
    #exec sudo -u -H R-ODAF -- bash
fi

#exec $@
#instead of sleeping, an HPC environment would be better suited to having
#an entrypoint that runs the pipeline; eaiser to do interactively, currently.
sleep infinity