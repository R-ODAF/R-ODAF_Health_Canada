# Metadata

This folder should contain at least the following file:  

metadata.txt  

The metadata file must be  a tab-delimited text file called "metadata.txt". Place this file in inputs/metadata/ .

Do not include spaces or special characters (i.e., anything except letters, numbers, dashes, or underscores) in the metadata file, as these will cause errors during analysis.

### Required columns
The following columns are required in the metadata file:
  - *sample_ID* (str) : The sample identifier column. No spaces or special characters allowed. Must be unique for each sample This identifier must match the fastq.gz file names (see below). 
  - *technical_control* (boolean): Allowed values are "T" or "F". Defines whether a sample was a technical control
  - *reference_rna* (boolean): Allowed values are "T" or "F". Defines whether a sample was an RNA control in a TempO-Seq experiment.  Set all to F for RNA-seq data.
  - *solvent_control* (boolean): Allowed values are "T" or "F". Defines whether a sample was a solvent control

### Optional columns
Any additional columns should be descriptive information about the samples. Suggested columns:

- *chemical*: The name of the treatment chemical
- *dose*: Treatment dose given
- *day*:  For time-series experiments
- *batch*: Batches are technical groups that could influence sequencing outcome. For example, if samples are sequenced in different library pools, include a column with pool number. If necessary, include multiple columns for batch variables (named appropriately)
- *cell_line*: if your experiment included multiple cell lines.

### Sample names
A column labeled *sample_ID* must be included in the metadata table, containing a unique identifier for each sample. This identifier must match the fastq.gz file names in the following way:

For single-end sequencing:  *{sample_ID}*.fastq.gz

For paired-end sequencing: *{sample_ID}*.R1.fastq.gz and *{sample_ID}*.R2.fastq.gz