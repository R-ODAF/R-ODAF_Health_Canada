# Metadata

This folder should contain at least the following file:  

metadata.txt  

The metadata file must be  a tab-delimited text file. Specify the file name in the config's *metadata_file* parameter.

Do not include spaces or special characters (i.e., anything except letters, numbers, dashes, or underscores) in the metadata file, as these will cause errors during analysis.

### Required columns
The following columns are required in the metadata file:
  - *sample_ID* (str) : The sample identifier column. No spaces or special characters allowed. Must be unique for each sample. This identifier must match the fastq.gz file names for [RNA-seq and TempO-seq experiments](#sample-names-in-rna-seq-and-tempo-seq-experiments). 
  - *technical_control* (boolean): Allowed values are "T" or "F". Defines whether a sample was a technical control (ex. lysis buffer only or )
  - *reference_rna* (boolean): Allowed values are "T" or "F". Defines whether a sample was an RNA control in a TempO-Seq experiment.  Set all to F for RNA-seq data.
  - *solvent_control* (boolean): Allowed values are "T" or "F". Defines whether a sample was an untreated control (commonly known as a solvent control or vehicle control in chemical exposure experiments)
Additional columns required for DRUG-Seq experiments:
  - *library_ID* (str): a library identifier column. Must match the input fastq file names. See [below](#sample-and-library-names-in-drug-seq-experiments)
  - *sample_barcode* (str): the per-sample unique barcode sequence. Used for assigning read counts to samples for the count table, and for demultiplexing reads into per-sample BAM files.

### Optional columns
Any additional columns should be descriptive information about the samples. Suggested columns:

- *chemical*: The name of the treatment chemical
- *dose*: Treatment dose given
- *day*:  For time-series experiments
- *hour*: for time-series experiments
- *batch*: Batches are technical groups that could influence sequencing outcome. For example, if samples are sequenced in different library pools, include a column with pool number. If necessary, include multiple columns for batch variables (named appropriately)
- *cell_line*: if your experiment included multiple cell lines.

### Sample names in RNA-seq and TempO-seq experiments
The *sample_ID* in each row must match the fastq.gz file names in the following way:

For single-end sequencing:  *{sample_ID}*.fastq.gz

For paired-end sequencing: *{sample_ID}*.R1.fastq.gz and *{sample_ID}*.R2.fastq.gz

### Sample and library names in DRUG-seq experiments
In DRUG-seq experiments, each pair of fastq files (produced by demultiplexing based on Illumina indices) equates to one library. Each library will contain multiple samples (which will be demultiplexed based on barcodes as part of this workflow).

The contents of the *library_ID* column must match the fastq files in the following way:

*{library_ID}*_R1.fastq.gz and *{library_ID}*_R2.fastq.gz

The contents of the *library_ID* column are not expected to be unique to each row, since each library contains multiple samples. You may only have a single library in your DRUG-Seq experiment.

The *sample_ID* in each row must be unique.