########################################################
## Build R-ODAF container for transcriptomic analysis ##
########################################################
# Mamba base image
FROM condaforge/mambaforge:23.3.1-1 as base
#continuumio/miniconda3 - difficult to install mamba, so switched away from this base 

# Required to avoid interactive prompts
ARG DEBIAN_FRONTEND=noninteractive
ARG PLATFORM=temposeq
ARG BUILD_CORES=8

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
      gosu nobody true  && \
      df -h && \
      apt-get clean && \
      df -h
      # verify that gosu works
      # free up disk space on runner

RUN mkdir -p /home/R-ODAF/R-ODAF_Health_Canada
WORKDIR "/home/R-ODAF/R-ODAF_Health_Canada"
COPY . .
RUN chown -R R-ODAF:R-ODAF /home/R-ODAF
USER R-ODAF
RUN mv inputs inputs.bak

RUN git clone https://github.com/EHSRB-BSRSE-Bioinformatics/test-data \
&& mv test-data/${PLATFORM}/* ./ \
&& wget https://github.com/EHSRB-BSRSE-Bioinformatics/unify_temposeq_manifests/raw/main/output_manifests/Human_S1500_1.2_standardized.csv \
&& df -h

# Build environments with Snakemake
RUN /bin/bash -c "snakemake --cores ${BUILD_CORES} --use-conda --conda-create-envs-only && conda clean -y -a && df -h"

# Install extra dependency for reports
RUN /bin/bash -c "conda run -p $(grep -rl "R-ODAF_reports" .snakemake/conda/*.yaml | sed s/\.yaml//) Rscript install.R"

# Print conda environment packages
RUN /bin/bash -c "\
  for env in .snakemake/conda/*; do \
    if [ -d \"$env\" ]; then \
      echo \"Contents of $env:\"; \
      source activate $env; \
      conda list; \
      conda deactivate; \
    fi; \
  done"

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
RUN /bin/bash -c "df -h && snakemake --cores 8 --use-conda"

# Clean up files
FROM tests as cleanup

# Clean up directories
RUN df -h &&\
      mkdir ./tests && \
      mv output \
      test-data \
      truth_checksums \
      wikipathways-20210810-gmt-Homo_sapiens.gmt \
      Human_S1500_1.2_standardized.csv \
      inputs \
      ./tests/
RUN /bin/bash -c "tree inputs.bak && \
                  mv inputs.bak inputs"
# Move test files
USER root
RUN rm -r ./tests

########################################################
## Build R-ODAF container for transcriptomic analysis ##
########################################################
# From tested container with updated code
FROM cleanup
ENTRYPOINT ["/home/R-ODAF/R-ODAF_Health_Canada/entrypoint.sh"]
