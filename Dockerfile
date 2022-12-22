########################################################
## Build R-ODAF container for transcriptomic analysis ##
########################################################
# Ubuntu base image
#FROM ubuntu:20.04
# Conda base image, test this instead
FROM continuumio/anaconda3

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
	libxt-dev \
	tree \
	vim

# Switch to analysis user
USER R-ODAF

# Install conda package manager to install tools (i.e., other dependencies)
#?ENV CONDA_DIR /home/R-ODAF/miniconda
#?RUN wget --quiet https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda.sh && \
#?	/bin/bash ~/miniconda.sh -b -p ${CONDA_DIR}
# Required for conda to be available on $PATH
#?ENV PATH=$CONDA_DIR/bin:$PATH
# Also install mamba for faster builds
RUN conda install -c conda-forge mamba

# Clone the R-ODAF repository
# Doing this at the build stage will mean that a given container is frozen for the version used here
# It should probably be recorded? Maybe use git hash?
RUN mkdir -p /home/R-ODAF/R-ODAF_Health_Canada
WORKDIR "/home/R-ODAF/R-ODAF_Health_Canada"
COPY . .
USER root
RUN chown -R R-ODAF:R-ODAF ./*
USER R-ODAF
RUN rm -r data \
&& rm -r config 

RUN git clone https://github.com/EHSRB-BSRSE-Bioinformatics/test-data \
&& mv test-data/temposeq/* ./ \
&& wget https://github.com/EHSRB-BSRSE-Bioinformatics/unify_temposeq_manifests/raw/main/output_manifests/Human_S1500_1.2_standardized.csv
# Load the conda environment
#RUN eval "$(conda shell.bash hook)"
RUN conda init
# Ensure that the reports environment has required R packages
SHELL ["conda", "run", "-n", "base", "/bin/bash", "-c"]

RUN conda activate base \
&& mamba create -c conda-forge -c bioconda -n snakemake snakemake \
&& conda activate snakemake \
&& snakemake --cores 8 --use-conda

# Install extra dependency
RUN conda activate R-ODAF_reports \
&& R -e "chooseCRANmirror(1, graphics=FALSE); remotes::install_github('bwlewis/crosstool')"

ENTRYPOINT ["/bin/bash", "-l", "-c"]
