########################################################
## Build R-ODAF container for transcriptomic analysis ##
########################################################
# Conda base image
FROM continuumio/miniconda3

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

#SHELL ["/bin/bash", "-c"]
RUN mkdir -p /home/R-ODAF/R-ODAF_Health_Canada
WORKDIR "/home/R-ODAF/R-ODAF_Health_Canada"
COPY . .
RUN conda install -n base -c conda-forge mamba
RUN mamba install -y -c conda-forge -c bioconda -n base snakemake
USER root
RUN rm -r data \
&& rm -r config
RUN chown -R R-ODAF:R-ODAF /home/R-ODAF
USER R-ODAF
RUN git clone https://github.com/EHSRB-BSRSE-Bioinformatics/test-data \
&& mv test-data/temposeq/* ./ \
&& wget https://github.com/EHSRB-BSRSE-Bioinformatics/unify_temposeq_manifests/raw/main/output_manifests/Human_S1500_1.2_standardized.csv

# Build environments with Snakemake
RUN snakemake --cores 32 --use-conda --conda-create-envs-only #--conda-prefix /home/R-ODAF/.conda/envs/
ARG FIND_HASH=$(grep -rl "R-ODAF_reports" .snakemake/conda/*.yaml | sed s/\.yaml//)
ENV ENV_HASH=${FIND_HASH}
# Install extra dependency for reports
SHELL ["conda", "run", "-p", "${ENV_HASH}", "/usr/bin/R", "install.R"]
#RUN R -e "chooseCRANmirror(1, graphics=FALSE); remotes::install_github('bwlewis/crosstool')"

RUN snakemake --cores 32

ENTRYPOINT ["/bin/bash", "-l", "-c"]
