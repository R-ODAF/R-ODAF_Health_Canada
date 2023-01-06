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
chown -R /home/R-ODAF/ R-ODAF:R-ODAF

# Drop privileges and execute next container command, or 'bash' if not specified.
if [[ $# -gt 0 ]]; then
    exec sudo -u -H R-ODAF -- "$@"
else
    exec sudo -u -H R-ODAF -- bash
fi
