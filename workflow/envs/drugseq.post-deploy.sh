#!/bin/bash
set -euxo pipefail

# Install FastReadCounter v1.1 or greater from https://github.com/DeplanckeLab/FastReadCounter
# As a single java jar
cd $CONDA_PREFIX/bin
wget https://github.com/DeplanckeLab/FastReadCounter/releases/download/v1.1/FastReadCounter-1.1.jar

echo "FastReadCounter 1.1 installed"