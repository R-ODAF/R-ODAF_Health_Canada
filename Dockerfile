########################################################
## Build R-ODAF container for transcriptomic analysis ##
########################################################
# Ubuntu base image
FROM ubuntu:20.04

# Required to avoid interactive prompts
ARG DEBIAN_FRONTEND=noninteractive

# Set up a user for running; -m creates home directory, -s sets default shell
RUN useradd -ms /bin/bash R-ODAF

# Install dependencies for running pipeline
RUN apt-get update && apt-get -y install \
        build-essential \
        wget \
        curl \
        git \
        libcairo2-dev \
        libxt-dev

# Install conda package manager to install tools (i.e., other dependencies)
RUN wget --quiet https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda.sh && \
        /bin/bash ~/miniconda.sh -b -p /opt/conda
# Required for conda to be available on $PATH
ENV PATH=/opt/conda/bin:$PATH
# Also install mamba for faster builds
RUN conda install -c conda-forge mamba

# Clone the R-ODAF repository
# Doing this at the build stage will mean that a given container is frozen for the version used here
# It should probably be recorded? Maybe use git hash?
# First, specify branch to use
USER R-ODAF
WORKDIR "/home/R-ODAF/"
COPY ./environment.yml ./environment.yml
COPY ./Rmd ./Rmd
COPY ./R-ODAF.Rproj ./R-ODAF.Rproj 
COPY ./scripts ./scripts
COPY ./workflow ./workflow
RUN git clone https://github.com/EHSRB-BSRSE-Bioinformatics/test-data \
&& mv test-data/temposeq/* ./ \
&& wget https://github.com/EHSRB-BSRSE-Bioinformatics/unify_temposeq_manifests/raw/main/output_manifests/Human_S1500_1.2_standardized.csv
#WORKDIR "/home/R-ODAF/R-ODAF_Health_Canada"
# Load the conda environment
RUN eval "$(conda shell.bash hook)"
# Run this if you don't already have an environment to use...
RUN mamba env create -f environment.yml
RUN conda init
SHELL ["conda", "run", "-n", "R-ODAF", "/bin/bash", "-c"]

#RUN conda activate R-ODAF

# Install extra dependencies
RUN R -e "chooseCRANmirror(1, graphics=FALSE); \
        remotes::install_github('bwlewis/crosstool'); \
        install.packages('cellWise')"

ENTRYPOINT ["/bin/bash"]