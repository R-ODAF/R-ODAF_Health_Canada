########################################################
## Build R-ODAF container for transcriptomic analysis ##
########################################################
# Mamba base image
FROM condaforge/mambaforge:23.3.1-1 as base
#continuumio/miniconda3 - difficult to install mamba, so switched away from this base 

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
        gosu \
	vim && \
      mamba install -y -c conda-forge -c bioconda -n base snakemake && \
      rm -rf /var/lib/apt/lists/* && \
      gosu nobody true 
      # last line is to verify that gosu works

RUN mkdir -p /home/R-ODAF/R-ODAF_Health_Canada
WORKDIR "/home/R-ODAF/R-ODAF_Health_Canada"
COPY . .
RUN chown -R R-ODAF:R-ODAF /home/R-ODAF
USER R-ODAF
RUN mv data data.bak && mv config config.bak

RUN git clone https://github.com/EHSRB-BSRSE-Bioinformatics/test-data \
&& mv test-data/temposeq/* ./ \
&& wget https://github.com/EHSRB-BSRSE-Bioinformatics/unify_temposeq_manifests/raw/main/output_manifests/Human_S1500_1.2_standardized.csv

# Build environments with Snakemake
RUN /bin/bash -c "snakemake --cores 32 --use-conda --conda-create-envs-only"

# Install extra dependency for reports
RUN /bin/bash -c "conda run -p $(grep -rl "R-ODAF_reports" .snakemake/conda/*.yaml | sed s/\.yaml//) Rscript install.R"

# Change ownership of files
USER root
RUN chown -R R-ODAF:R-ODAF /home/R-ODAF
USER R-ODAF

########################################################
## Build R-ODAF container for transcriptomic analysis ##
########################################################
# R-ODAF base image
# Run tests using base dependencies, but updated code
# Note that if you alter any dependencies, the base container must be rebuilt
FROM base as tests

# Run tests
RUN /bin/bash -c "snakemake --cores 32 --use-conda"

# Clean up files
FROM tests as cleanup

# Clean up directories
RUN mkdir ./tests && \
    mv analysis \
       test-data \
       truth_checksums \
       wikipathways-20210810-gmt-Homo_sapiens.gmt \
       Human_S1500_1.2_standardized.csv \
       logs \
       data \
       config \
       ./tests/
RUN /bin/bash -c "tree data.bak && tree config.bak && \
                  mv data.bak data && mv config.bak config"
# Move test files
USER root
RUN mv ./tests /opt/tests

########################################################
## Build R-ODAF container for transcriptomic analysis ##
########################################################
# From tested container with updated code
FROM cleanup
ENTRYPOINT ["/home/R-ODAF/R-ODAF_Health_Canada/entrypoint.sh"]
