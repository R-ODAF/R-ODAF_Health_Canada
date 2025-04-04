########################################################
## Build R-ODAF container for transcriptomic analysis ##
########################################################
# Mamba base image
# snakemake image is based on mambaorg/micromamba
# mambaorg/micromamba in turn is based on debian:bookworm-slim
FROM snakemake/snakemake:stable
# v
# v8.4.8


# Required to avoid interactive prompts
#ARG DEBIAN_FRONTEND=noninteractive
ARG PLATFORM=temposeq
ARG BUILD_CORES=8
RUN useradd -ms /bin/bash R-ODAF
# Install runtime dependencies that are necessary for your application
# Why /var/lib/dpkg/lock-frontend doesn't exist... I have no idea
RUN mkdir -p /var/lib/dpkg && touch /var/lib/dpkg/lock-frontend \
      &&apt-get update \
      && apt-get install --no-install-recommends gosu  \
      && gosu nobody true \
      && rm -rf /var/lib/apt/lists/* \
      && apt-get clean \
      && df -h

RUN mkdir -p /home/R-ODAF/R-ODAF_Health_Canada
WORKDIR "/home/R-ODAF/R-ODAF_Health_Canada"
COPY . .
RUN mv inputs inputs.bak \
      && wget https://github.com/EHSRB-BSRSE-Bioinformatics/test-data/archive/master.tar.gz \
      && tar -zxvf master.tar.gz \
      && mv -f test-data-main/${PLATFORM}/* ./ \
      && wget https://github.com/EHSRB-BSRSE-Bioinformatics/unify_temposeq_manifests/raw/main/output_manifests/Human_S1500_1.2_standardized.csv \
      && rm master.tar.gz

RUN chown -R R-ODAF:R-ODAF /home/R-ODAF && df -h
USER R-ODAF

# Build environments with Snakemake
RUN snakemake --cores ${BUILD_CORES} --software-deployment-method conda --conda-create-envs-only 

# Install extra dependency for reports
RUN conda run -p $(grep -rl "R-ODAF_reports" .snakemake/conda/*.yaml | sed s/\.yaml//) Rscript install.R 

RUN snakemake --cores ${BUILD_CORES} --software-deployment-method conda \
      && df -h \
      && mkdir ./tests \
      && mv output test-data-main truth_checksums wikipathways-20210810-gmt-Homo_sapiens.gmt Human_S1500_1.2_standardized.csv inputs ./tests/ \
      && mv inputs.bak inputs \
      && rm -rf ./tests

USER root
# Set the entrypoint for the container
ENTRYPOINT ["/home/R-ODAF/R-ODAF_Health_Canada/entrypoint.sh"]
