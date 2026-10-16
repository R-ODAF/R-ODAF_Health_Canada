# Reference files

This folder can contain the genome reference files (fasta, .gtf, and, if necessary, biospyder manifest files or ercc fasta and gtf files). 

Reference files can also be stored elsewhere on your computer. We recommend doing this if you will be doing multiple projects, so that reference files are not replicated on your file system.

The path to your reference files must be given in the config file (in **genomedir**, **biospyder_dbs**, and **erccdir** parameters).
The STAR index will be made in the same location as the reference files.

The names of your reference files must be given in the config file (**genome_filename**, **annotation_filename**, **biospyder_manifest_file**, )

## RNA-seq and DRUG-seq: fasta and gtf files
Provide a genome fasta (for example, Homo_sapiens.GRCh38.dna.primary_assembly.fa available from https://www.gencodegenes.org/human) and an associated annotation file in gtf format (ex. Homo_sapiens.GRCh38.116.gtf.gz available from https://ftp.ensembl.org/pub/release-116/gtf/homo_sapiens/). 

NOTE: Some annotation files (gtf files) include a version number appended to the end of EnsemblIDs (ex. "ENSG00000290825.1"). This is not compatible with the code in this workflow. Please use an annotation file without version numbers.

The names of these reference files must be given in the config file (**genome_filename** and **annotation_filename** parameters).

## DRUG-Seq (optional) ERCC files
If you included ERCC spike-ins during library building, you can set the config parameter **include_ercc** to TRUE to check the percent of mapped reads matching to ERCC transcripts. You must supply ERCC fasta and gtf files, which will be concatenated with the genome fasta and gtf files for your organism.

## TempO-Seq manifest, fasta, and gtf files
TempO-seq probe sequences and annotations are derived from the manifest files provided by BioSpyder. We have standardized the format of the most commonly used TempO-seq kits, with manifests available here: https://github.com/EHSRB-BSRSE-Bioinformatics/unify_temposeq_manifests/.

You can produce fasta and gtf files from a BioSpyder manifest using the biospyder_manifest_to_reffiles.py script in the scripts directory.
