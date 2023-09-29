########################################################
## Build R-ODAF container for transcriptomic analysis ##
########################################################
# Conda base image
FROM continuumio/miniconda3 as base 
#:4.12.0 as base

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

RUN mkdir -p /home/R-ODAF/R-ODAF_Health_Canada
WORKDIR "/home/R-ODAF/R-ODAF_Health_Canada"
COPY . .
RUN conda install 'mamba<=1.4.5' -n base -c conda-forge
RUN conda update mamba -n base -c conda-forge
#RUN conda install -n base -c conda-forge mamba
RUN mamba install -y -c conda-forge -c bioconda -n base snakemake
RUN chown -R R-ODAF:R-ODAF /home/R-ODAF
USER R-ODAF
RUN mv data data.bak && mv config config.bak

RUN git clone https://github.com/EHSRB-BSRSE-Bioinformatics/test-data \
&& mv test-data/temposeq/* ./ \
&& wget https://github.com/EHSRB-BSRSE-Bioinformatics/unify_temposeq_manifests/raw/main/output_manifests/Human_S1500_1.2_standardized.csv

# Build environments with Snakemake
RUN snakemake --cores 32 --use-conda --conda-create-envs-only

# Install extra dependency for reports
RUN conda run -p $(grep -rl "R-ODAF_reports" .snakemake/conda/*.yaml | sed s/\.yaml//) Rscript install.R

########################################################
## Build R-ODAF container for transcriptomic analysis ##
########################################################
# R-ODAF base image
# Run tests using base dependencies, but updated code
# Note that if you alter any dependencies, the base container must be rebuilt
FROM base as tests

# Update code if it has changed in build context
COPY . . 
USER root
RUN chown -R R-ODAF:R-ODAF /home/R-ODAF
USER R-ODAF

RUN /bin/bash -c "snakemake --cores 32 --use-conda"

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
RUN sleep 1 && mv -vf "data.bak" "data" && \
    mv -vf "config.bak" "config"

########################################################
## Build R-ODAF container for transcriptomic analysis ##
########################################################
# From tested container with updated code
FROM cleanup

USER root

# Move files
RUN mv tests /opt/tests

# Install gosu
RUN set -eux; \
	apt-get update; \
	apt-get install -y gosu; \
	rm -rf /var/lib/apt/lists/*; \
# verify that the binary works
	gosu nobody true

COPY entrypoint.sh .
RUN ["chmod", "+x", "entrypoint.sh"]
ENTRYPOINT ["/home/R-ODAF/R-ODAF_Health_Canada/entrypoint.sh"]
#CMD ["/bin/bash" "-l"]
