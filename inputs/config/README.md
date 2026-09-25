# Config file

Information from the config file is used throughout this pipeline, both by the Snakemake workflows and in the R scripts that run analyses and generate reports. 

A default file (inputs/config/configs.default.yaml) is provided that contains examples of all the necessary parameters. To copy the default file into a config.yaml file that you can populate with the correct information for your project, run `install.sh` .



## Parameters

## Metadata and config files
Provide the names of the metadata and contrast files to use in this analysis. These files must be located in inputs/metadata/ and inputs/contrasts, respectively. 
This allows the user flexilibity to, for example, run the differential expression step with one set of contrasts, then edit the config file to specify a different contrasts file and re-run module 3 of the workflow to do an additional analysis.

Note that the contrasts file used for each differential expression analysis is copied into output/analysis/{analysis_dir_name}/Pipeline_record.

### Batch variable

### Manifest file (for Tempo-seq experiments)
Biospyder provides manifest files for each TempO-seq kit, but they do not follow a standard format and are often missing information such as EntrezID or EnsemblID. To address this, we have created standardized TempO-Seq manifests for frequently used kits, available for download from https://github.com/EHSRB-BSRSE-Bioinformatics/unify_temposeq_manifests/tree/main/output_manifests


### Genome

You must provide a reference sequence for alignment of the raw reads. For TempO-Seq experiments, these are provided by BioSpyder; for RNA-seq experiments they can be obtained through your favorite database. This is currently beyond the scope of this guide. Ensure you have some type of annotation file (GTF format) available as well, to dictate which sequences in the FASTA file correspond to which genes or probes. These should be used to create a STAR index within the directory where you store your reference genome in FASTA format.


### Additional DESeq2 parameters

group_facet:
Setting this to a column name will cause separate DESeq2 analyses to be done on each group in that column. For example, if you have a column called "chemicals" containing A, B, and C, and you set group_facet: "chemicals", separate DESeq2 object will be made for A, B, and C.


* **Strict and Lenient contrasts** During an analysis in which experiments are faceted with `group_facet` and filtered with `deseq_filter`, contrasts are filtered to match the entries in `deseq_filter`. If neither `strict_contrasts` or `lenient_contrasts` is true, then only the experimental element of a contrast is tested for membership in `deseq_filter`. If `strict_contrasts` is `TRUE`, both the experimental and control elements of a contrast must be in `deseq_filter` for that contrast to be examined. If `lenient_contrasts` is `TRUE`, either one of the experimental or control elements is enough for a contrast to be included.