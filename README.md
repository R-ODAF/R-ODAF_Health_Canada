# Omics Data Analysis Frameworks for Regulatory application (R-ODAF)

## What is the R-ODAF?

It is an analysis framework geared toward the use of 'omics data in toxicology experiments. It aims to use commonly available open-source tools for analysis of transcriptomic data, specifically using DESeq2 for the determination of differentially expressed genes (DEGs).

This R-ODAF repository, which was forked from the main branch, is under development by researchers at Health Canada and has been used primarily in the context of regulatory toxicology. The intended use is for processing high-throughput sequencing data from FASTQ files to DEG lists, and to generate exploratory analyses for use in toxicology research.

## What does this repository do?

This repository can be used to generate:

* study-wide HTML reports on quality of the HTS data using MultiQC
* study-wide HTML reports on the quality of alignments and biosamples within the study using various metrics
* lists of DEGs
  * DEG lists can be sliced several different ways depending on the application
  * Downstream application of BMD modeling is an option
  * Tools such as IPA can be used to do further downstream analysis
  * R-ODAF filtering criteria as set out in publications is applied to DEG lists
* exploratory data visualization for DEGs, inspired by the regionReport DESeq2 reports
* interactive tables and plots related to analysis of gene expression

## What does this repository NOT do?

This repository provides no warranty or guarantees on the quality of your analysis. It is used actively in research projects and is not intended to provide any type of authoritative assessment of your data.

There may be other tools that are more applicable to your goals. Indeed, if you are not undertaking toxicogenomics experiments, then that is likely the case. 

## How can I use this repository?

In its current state of development, this repository should be used by individuals (or teams) experienced with running R in a Linux-based environment. It is necessary to understand scripting and how to process data. We encourage users to learn if they are interested in pursuing this type of data analysis; to that end, we provide a users guide below detailing the steps required to process data, and best practices for using this analysis pipeline.

# Getting Started

## Cloning the repository

You should start by cloning the repository locally, to a computer that can handle analysis of your data. Preferably, one where your raw data is locally acessible.

You can clone the repository using this command: 

`git clone https://github.com/R-ODAF/R-ODAF_Health_Canada.git`

## Working in RStudio

Typically, when undertaking a new analysis, it is a good idea to have a separate folder for each project. This makes it easier to manage, and also makes your code more portable. Working in RStudio facilitates these tasks.

Therefore, it is a good idea to open the cloned repository in RStudio as a project (there is already a project file by default in this repository, so RStudio should recognize it). This tells RStudio where to look for files, which is better than setting a working directory.

## Prerequisites for running analyses

There are many R scripts in this repository, each with different purposes, as well as a portion of the pipeline that is run using [Snakemake](https://snakemake.readthedocs.io/en/stable/). I suggest reading over that documentation if you are not already familiar with it, as it is beyond the scope of this documentation.

### List of dependencies

The following software is required to pre-process high-throughput sequencing data:

* samtools
* STAR
* MultiQC
* fastp

It is a good strategy to manage these particular software dependencies, required on your path, by setting up a Conda environment with them installed. 

There are also numerous dependencies in R. This list might be complete: 

AnnotationHub
BiocParallel
biomaRt
Cairo
cellWise
cluster
clusterProfiler
ComplexHeatmap
conflicted
crosstalk
crosstool
data.table
dendextend
dendsort
DESeq2
doParallel
DT
edgeR
enrichplot
fields
foreach
GGally
ggplot2
ggridges
gtools
heatmap3
kableExtra
knitr
lattice
limma
magrittr
manipulateWidget
openxlsx
org.(YourFavoriteGenomeHere).eg.db
OUTRIDER
pheatmap
plotly
RColorBrewer
RMariaDB
rrcov
rWikiPathways
sessioninfo
sfsmisc
tidytext
tidyverse
UpSetR
viridis
vsn
vtree
yaml

# Running the data analysis

## "Installing" the pipeline

In the top level directory of the cloned repository, you can run the small utility install.sh (i.e., in a bash terminal, you can run `bash ./install.sh`). This will copy the default configuration file

The “install.sh” script will create a config.yaml file, which you will need to customize (see below).

## Setting up your configuration

### Populate config.yaml

Once you have `config.yaml`, this file will dictate almost everything you need to customize about your experiment (with the exception of the input data itself). Go through and fill out the necessary parameters. It is used by the Snakemake workflows, as well as the render_ R scripts to generate reports and run analysis code.

### Genome parameters

You must provide a reference sequence for alignment of the raw reads. For TempO-Seq experiments, these are provided by BioSpyder; for RNA-seq experiments they can be obtained through your favorite database. This is currently beyond the scope of this guide. Ensure you have some type of annotation file (GTF format) available as well, to dictate which sequences in the FASTA file correspond to which genes or probes. These should be used to create a STAR index within the directory where you store your reference genome in FASTA format.

### Metadata

You must provide a metadata file (sometimes also known as 'colData', i.e., data about the columns of your count matrix; sample metadata; or sample information) describing the experimental parameters.

### Contrasts

You must provide a simple text file to describe which contrasts of interest should be tested with the `results()` function of DESeq2.

## Pre-processing

The first step is running pre-processing of FASTQ files using the Snakemake pipeline. It’s the workflow/Snakemake.temposeq file that is run for the pipeline.

Run this using a command like this:

`snakemake --cores 40 --snakefile ./workflow/Snakefile.temposeq`

## Study-wide Quality Control

The `./Rmd/Sample_QC.Rmd` file is used to generate an HTML report on the quality of your alignments.

Rather than knitting this within RStudio, it is preferable to call the `render()` function of `RMarkdown` using `./scripts/render_studywide_QC_report.R`. This will then use `config.yaml` to populate any custom `params` in `RMarkdown` as required.

## Running DESeq2 and Generating Exploratory Analysis

Once you have run the study-wide QC, it will output a new metadata file (`metadata.QC_applied.txt`) which will only contain samples passing the QC filtering. You can see which samples were removed in `./analysis/QC/details`.

Again, rather than using the `./Rmd/DESeq2_report.rnaseq.Rmd` file interactively or knitting it in RStudio, it is preferable to wrap the rendering of reports using `./scripts/render_DESeq2_report.R`, which will also handle splitting things up for large experiments. For example, if you have a study with 10 chemicals, each of which produces gene expression profiles that are vastly different, then you likely want to analyze these chemicals separately. You might want a report showing all chemicals. You might want a report showing only one chemical. These options are both possible by specifying parameters in the `config.yaml` file.

## Questions?

I'm not surprised if you have questions. This is very much a draft of the workflow, and is incomplete except for serving highly experienced users. Please reach out if you have questions or issues, and feel free to submit an [issue](https://github.com/R-ODAF/R-ODAF_Health_Canada/issues) or PR with any suggestions you might have.
