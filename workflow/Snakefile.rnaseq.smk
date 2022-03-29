import re
import pandas as pd
from glob import glob
import pathlib
from snakemake.utils import validate

# load and validate config stuff
configfile: "config/config.yaml"

# set up directories
main_dir = Path(config["prefix"])

genome_dir = Path(config["genomedir"])

data_dir = main_dir / "data"
raw_dir = data_dir / "raw" # should already exist
metadata_dir = data_dir / "metadata" # should already exist
metadata_file = metadata_dir / "metadata.txt" # should aready exist
#TODO set up a check to make sure this stuff exists


SAMPLES = pd.read_table(metadata_file)['sample_id'].tolist()

print("samples: " + str(SAMPLES))

# set up output dirs

processed_dir = data_dir / "processed"
processed_dir.mkdir(parents=True, exist_ok=True)

trim_dir = processed_dir / "trimmed"
align_dir = processed_dir / "aligned"
quant_dir = processed_dir / "quant"
trim_dir.mkdir(parents=True, exist_ok=True)
align_dir.mkdir(parents=True, exist_ok=True)
quant_dir.mkdir(parents=True, exist_ok=True)

analysis_dir = main_dir / "analysis"
analysis_dir.mkdir(parents=True, exist_ok=True)

qc_dir = analysis_dir / "QC"
qc_dir.mkdir(parents=True, exist_ok=True)

log_dir = main_dir / "logs"
log_dir.mkdir(parents=True, exist_ok=True)

print("using genome: " + str(genome_dir))

shell("cp -rf workflow {log_dir}")
shell("cp -rf config {log_dir}")

# run whole pipeline

rule all:
    input:
        expand(str(trim_dir / "{sample}.fastq.gz"), sample=SAMPLES),
        qc_dir / "multiqc_report.html",
        processed_dir / "genes.data.tsv",
        processed_dir / "isoforms.data.tsv" 
    log: log_dir / "all.log"


##################################
### Trimming raw reads : Fastp ###
##################################

################################
#INFORMATION ON TRIMMING PROCESS
#--cut_front --cut_front_window_size 1 --cut_front_mean_quality 3 == Trimmomatic  "LEADING:3"
#--cut_tail --cut_tail_window_size 1 --cut_tail_mean_quality 3 == Trimmomatic  "TRAILING:3"
#--cut_right --cut_right_window_size 4 --cut_right_mean_quality 15 == Trimmomatic  "SLIDINGWINDOW:4:15"

#### If additional trimming is needed (see multiQCreport): 
# add to R1: --trim_front1 {amount_bases} and/or --trim_tail1 {amount_bases}
# add to R2: --trim_front2 {amount_bases} and/or --trim_tail2 {amount_bases}


rule fastp_all:
    input:
        expand(str(trim_dir / "{sample}.fastq.gz"), sample=SAMPLES),
    log: log_dir / "bbduk_all.log"

if config["mode"] == "pe":
    rule fastp:
        input:
            R1=ancient(str(raw_dir / "{sample}.R1.fastq.gz")),
            R2=ancient(str(raw_dir / "{sample}.R2.fastq.gz")),
        output:
            R1 = trim_dir / "{sample}.R1.fastq.gz",
            R2 = trim_dir / "{sample}.R2.fastq.gz",
            json = trim_dir / "{sample}_fastp.json",
            html = trim_dir / "{sample}_fastp.html",
        log:
            filter_stats = log_dir / "bbduk.{sample}.log"
        params:
            front_size = 1,
            front_quality = 3,
            tail_size = 1,
            tail_quality = 3,
            right_size = 4,
            right_quality = 15,
            length_required = 36,
            threads = config["threads"],
        shell:
            '''
            fastp \
            --in1 {input.R1} \
            --in2 {input.R2} \
            --out1 {output.R1} \
            --out2 {output.R2} \
            --json {output.json} \
            --html {output.html} \
            --cut_front \
            --cut_front_window_size {params.front_size} \
            --cut_front_mean_quality {params.front_quality} \
            --cut_tail \
            --cut_tail_window_size {params.tail_size} \
            --cut_tail_mean_quality {params.tail_quality} \
            --cut_right \
            --cut_right_window_size {params.right_size} \
            --cut_right_mean_quality {params.right_quality} \
            --length_required {params.length_required} \
	    	--thread {params.threads} \
            '''

if config["mode"] == "se":
    rule fastp:
        input:
            R1=ancient(str(raw_dir / "{sample}.fastq.gz")),
        output:
            R1 = trim_dir / "{sample}.fastq.gz",
            json = trim_dir / "{sample}_fastp.json",
            html = trim_dir / "{sample}_fastp.html",
        log:
            filter_stats = log_dir / "bbduk.{sample}.log"
        params:
            front_size = 1,
            front_quality = 3,
            tail_size = 1,
            tail_quality = 3,
            right_size = 4,
            right_quality = 15,
            length_required = 36,
            threads = config["threads"],
        shell:
            '''
            fastp \
            --in1 {input.R1} \
            --out1 {output.R1} \
            --json {output.json} \
            --html {output.html} \
            --cut_front \
            --cut_front_window_size {params.front_size} \
            --cut_front_mean_quality {params.front_quality} \
            --cut_tail \
            --cut_tail_window_size {params.tail_size} \
            --cut_tail_mean_quality {params.tail_quality} \
            --cut_right \
            --cut_right_window_size {params.right_size} \
            --cut_right_mean_quality {params.right_quality} \
            --length_required {params.length_required} \
	    	--thread {params.threads} \
            '''


####################################################
### Alignment of reads (paired end & single end) ###
####################################################

rule STAR_make_index:
    input:
        genome = ancient(genome_dir / config["genome_filename"])
    params:
        index_dir = genome_dir / "STAR_index",
        annotations = genome_dir / config["annotation_filename"],
        overhang = 100,
        threads = config["threads"],
        suffix_array_sparsity = 2, # bigger for smaller (RAM), slower indexing
        genomeChrBinNbits = 12, # might need to mess with this for different genomes
    output:
        directory(genome_dir / "STAR_index"),
    shell:
        '''
		STAR \
        --runMode genomeGenerate \
		--genomeDir {params.index_dir} \
		--genomeFastaFiles {input.genome} \
		--sjdbGTFfile {params.annotations} \
		--sjdbOverhang {params.overhang} \
		--runThreadN {params.threads} \
		--genomeSAsparseD {params.suffix_array_sparsity} \
		--genomeChrBinNbits {params.genomeChrBinNbits}
        '''

rule STAR_all:
    input:
        expand(str(align_dir / "{sample}.Aligned.toTranscriptome.out.bam"), sample=SAMPLES)

if config["mode"] == "pe":
    rule STAR:
        input:
            R1 = trim_dir / "{sample}.R1.fastq.gz",
            R2 = trim_dir / "{sample}.R2.fastq.gz",
            index = genome_dir / "STAR_index",
        output:
            sortedByCoord = align_dir / "{sample}.Aligned.sortedByCoord.out.bam",
            toTranscriptome = align_dir / "{sample}.Aligned.toTranscriptome.out.bam"
        params:
            annotations = genome_dir / config["annotation_filename"],
            threads = config["threads"],
            bam_prefix = lambda wildcards : align_dir / "{}.".format(wildcards.sample),
        resources:
            load=100
        shell:
            '''
		    STAR \
			--runThreadN {params.threads} \
			--genomeDir {input.index} \
			--readFilesIn {input.R1} {input.R2} \
			--quantMode TranscriptomeSAM \
			--readFilesCommand zcat \
			--outFileNamePrefix {params.bam_prefix} \
            --outSAMtype BAM SortedByCoordinate
            '''

if config["mode"] == "se":
    rule STAR:
        input:
            R1 = trim_dir / "{sample}.fastq.gz",
            index = genome_dir / "STAR_index",
        output:
            sortedByCoord = align_dir / "{sample}.Aligned.sortedByCoord.out.bam",
            toTranscriptome = align_dir / "{sample}.Aligned.toTranscriptome.out.bam"
        params:
            annotations = genome_dir / config["annotation_filename"],
            threads = config["threads"],
            bam_prefix = lambda wildcards : align_dir / "{}.".format(wildcards.sample),
        resources:
            load=100
        shell:
            '''
            STAR \
                --runThreadN {params.threads} \
			    --genomeDir {input.index} \
                --readFilesIn {input.R1} \
                --quantMode TranscriptomeSAM \
                --readFilesCommand zcat \
                --outFileNamePrefix {params.bam_prefix} \
                --outSAMtype BAM SortedByCoordinate
            '''


#######################
# QUANTIFICATION RSEM #
#######################

rule RSEM_make_index:
    input:
        genome = genome_dir / config["genome_filename"]
    params:
        annotations = genome_dir / config["annotation_filename"],
        genome_name = config["genome_name"],
    output:
        directory(genome_dir / "RSEM_index"),
    shell:
        '''
        mkdir RSEM_index
        rsem-prepare-reference --gtf {params.annotations} {input.genome} {output}/{params.genome_name}
        '''

if config["mode"] == "pe":
    rule RSEM:
        input:
            bam = align_dir / "{sample}.Aligned.toTranscriptome.out.bam",
            index = genome_dir / "RSEM_index" 
        output:
            isoforms = quant_dir / "{sample}.isoforms.results",
            genes = quant_dir / "{sample}.genes.results",
        params:
            threads = config["threads"],
            output_prefix =  lambda wildcards : quant_dir / "{}".format(wildcards.sample),
            index_name = config["genome_name"],
        shell:
            '''
            rsem-calculate-expression \
            -p {params.threads} \
            --paired-end \
            --bam {input.bam} \
            --no-bam-output \
            {input.index}/{params.index_name} \
            {params.output_prefix}
            '''

if config["mode"] == "se":
    rule RSEM:
        input:
            bam = align_dir / "{sample}.Aligned.toTranscriptome.out.bam",
            index = genome_dir / "RSEM_index" 
        output:
            isoforms = quant_dir / "{sample}.isoforms.results",
            genes = quant_dir / "{sample}.genes.results",
        params:
            threads = config["threads"],
            output_prefix =  lambda wildcards : quant_dir / "{}".format(wildcards.sample),
            index_name = config["genome_name"],
        shell:
            '''
            rsem-calculate-expression \
            -p {params.threads} \
            --bam {input.bam} \
            --no-bam-output \
            {input.index}/{params.index_name} \
            {params.output_prefix}
            '''

rule counts_matrix:
    input:
        genes = expand(str(quant_dir / "{sample}.genes.results"), sample=SAMPLES),
        isoforms = expand(str(quant_dir / "{sample}.isoforms.results"), sample=SAMPLES)
    output:
        genes = processed_dir / "genes.data.tsv",
        isoforms = processed_dir / "isoforms.data.tsv"
    shell:
        '''
        rsem-generate-data-matrix {input.genes} > {output.genes}
        sed -i 's/\.genes.results//g' {output.genes}
        sed -i 's|{quant_dir}/||g' {output.genes}
        sed -i 's/"//g' {output.genes}
        rsem-generate-data-matrix {input.isoforms} > {output.isoforms}
        sed -i 's/\.isoforms.results//g' {output.isoforms}
        sed -i 's|{quant_dir}/||g' {output.isoforms}
        sed -i 's/"//g' {output.isoforms}
        '''




###############
### MultiQC ###
###############

rule multiqc:
    message: "running multiqc"
    input:
        expand(str(trim_dir / "{sample}_fastp.json"), sample=SAMPLES),
        expand(str(align_dir / "{sample}.Aligned.toTranscriptome.out.bam"), sample=SAMPLES),
        expand(str(quant_dir / "{sample}.genes.results"), sample=SAMPLES),
        expand(str(quant_dir / "{sample}.isoforms.results"), sample=SAMPLES)
    output:
        qc_dir / "multiqc_report.html"
    log: log_dir / "multiqc.log"
    run:
        shell("multiqc --cl_config \"extra_fn_clean_exts: {{ '_fastp.json' }}\" -fz {raw_dir} {trim_dir} {align_dir} {quant_dir}")
        shell("mv multiqc_report.html multiqc_data.zip {qc_dir}")
