# Config file

### Batch variable parameter

### Manifest file parameter (for Tempo-seq experiments)



### Genome parameters

You must provide a reference sequence for alignment of the raw reads. For TempO-Seq experiments, these are provided by BioSpyder; for RNA-seq experiments they can be obtained through your favorite database. This is currently beyond the scope of this guide. Ensure you have some type of annotation file (GTF format) available as well, to dictate which sequences in the FASTA file correspond to which genes or probes. These should be used to create a STAR index within the directory where you store your reference genome in FASTA format.


### Additional DESeq2 parameters

group_facet:
Setting this to a column name will cause separate DESeq2 analyses to be done on each group in that column. For example, if you have a column called "chemicals" containing A, B, and C, and you set group_facet: "chemicals", separate DESeq2 object will be made for A, B, and C.


* **Strict and Lenient contrasts** During an analysis in which experiments are faceted with `group_facet` and filtered with `group_filter`, contrasts are filtered to match the entries in `group_filter`. If neither `strict_contrasts` or `lenient_contrasts` is true, then only the experimental element of a contrast is tested for membership in `group_filter`. If `strict_contrasts` is `TRUE`, both the experimental and control elements of a contrast must be in `group_filter` for that contrast to be examined. If `lenient_contrasts` is `TRUE`, either one of the experimental or control elements is enough for a contrast to be included.