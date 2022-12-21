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

# Switch to analysis user
USER R-ODAF

# Install conda package manager to install tools (i.e., other dependencies)
ENV CONDA_DIR /home/R-ODAF/miniconda
RUN wget --quiet https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda.sh && \
	/bin/bash ~/miniconda.sh -b -p ${CONDA_DIR}
# Required for conda to be available on $PATH
ENV PATH=$CONDA_DIR/bin:$PATH
# Also install mamba for faster builds
RUN conda install -c conda-forge mamba

# Clone the R-ODAF repository
# Doing this at the build stage will mean that a given container is frozen for the version used here
# It should probably be recorded? Maybe use git hash?

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
# Set working directory to home
# First, specify branch to use
ARG BRANCH="main"
# Set working directory to home
RUN ls -alht
WORKDIR "/home/R-ODAF/"
RUN ls -alht
RUN git clone https://github.com/R-ODAF/R-ODAF_Health_Canada.git \
	&& cd R-ODAF_Health_Canada \
	&& git checkout ${BRANCH} \
	&& git clone https://github.com/EHSRB-BSRSE-Bioinformatics/test-data \
	&& rm -r data \ 
	&& rm -r config \
	&& mv test-data/temposeq/* ./ \
	&& wget https://github.com/EHSRB-BSRSE-Bioinformatics/unify_temposeq_manifests/raw/main/output_manifests/Human_S1500_1.2_standardized.csv
RUN ls -alht
WORKDIR "/home/R-ODAF/R-ODAF_Health_Canada"
# Load the conda environment
RUN eval "$(conda shell.bash hook)"
# Run this if you don't already have an environment to use...
RUN echo ${HOME}/miniconda/etc/profile.d/conda.sh >> ~/.bashrc
RUN mamba env create -f environment.yml

RUN echo "conda activate R-ODAF" >> ~/.bashrc
SHELL ["conda", "run", "-n", "R-ODAF", "/bin/bash", "-c"]
RUN echo "Checking if python is installed..."
RUN python -c "print('hello')" && conda info && conda list
RUN ls -alht
RUN echo "Checking if STAR is installed..."
RUN STAR --version
RUN R -e "chooseCRANmirror(1, graphics=FALSE); remotes::install_github('bwlewis/crosstool'); install.packages('cellWise')"

RUN snakemake --cores 8
RUN echo "bash ${HOME}/miniconda/etc/profile.d/conda.sh" >> ~/.bashrc
RUN conda init bash
CMD conda activate R-ODAF

RUN Rscript scripts/render_studywide_QC_report.R
RUN Rscript scripts/run_DESeq2.R
#RUN Rscript scripts/render_DESeq2_report.parallel.R

