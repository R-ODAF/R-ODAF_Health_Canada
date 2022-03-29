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

temposeq_script = main_dir / "scripts/temposeq/TempO-SeqR_v3.1.R"
temposeq_quantification_script = main_dir / "scripts/temposeq/temposeq_quantification.R"

sample_id_col = config["sample_id"]

SAMPLES = pd.read_table(metadata_file)[sample_id_col].tolist()

print("samples: " + str(SAMPLES))

# set up output dirs

processed_dir = data_dir / "processed"
processed_dir.mkdir(parents=True, exist_ok=True)

trim_dir = processed_dir / "trimmed"
align_dir = processed_dir / "aligned"
trim_dir.mkdir(parents=True, exist_ok=True)
align_dir.mkdir(parents=True, exist_ok=True)

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
        qc_dir / "MultiQC_Report.html",
        processed_dir / "count_table.csv",
        processed_dir / "mapped_unmapped.csv"
    log: log_dir / "temposeq.log"


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


rule fastp:
    input:
        R1=ancient(str(raw_dir / "{sample}.fastq.gz")),
    output:
        R1 = pipe(trim_dir / "{sample}.fastq.gz"),
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
        length_required = 50,
        trim_tail1 = 1,
        threads = config["threads"],
    shell:
        '''
        fastp \
        --in1 {input.R1} \
        --out1 {output.R1} \
        --json {output.json} \
        --html {output.html} \
        --cut_front \
        --disable_adapter_trimming \
        --cut_front_window_size {params.front_size} \
        --cut_front_mean_quality {params.front_quality} \
        --cut_tail \
        --cut_tail_window_size {params.tail_size} \
        --cut_tail_mean_quality {params.tail_quality} \
        --cut_right \
        --cut_right_window_size {params.right_size} \
        --cut_right_mean_quality {params.right_quality} \
        --length_required {params.length_required} \
        --trim_tail1 {params.trim_tail1} \
        --thread {params.threads} \
        '''

####################################################
### Alignment of reads (paired end & single end) ###
####################################################


rule STAR_all:
    input:
        expand(str(align_dir / "{sample}.Aligned.toTranscriptome.out.bam"), sample=SAMPLES)



rule STAR:
    input:
        R1 = trim_dir / "{sample}.fastq.gz",
        index = genome_dir,
    output:
        sortedByCoord = align_dir / "{sample}.Aligned.sortedByCoord.out.bam",
        toTranscriptome = align_dir / "{sample}.Aligned.toTranscriptome.out.bam"
    params:
        annotations = genome_dir / config["annotation_filename"],
        threads = config["threads"],
        bam_prefix = lambda wildcards : align_dir / "{}.".format(wildcards.sample),
    shell:
        '''
        STAR \
            --alignEndsType EndToEnd \
            --genomeLoad LoadAndKeep \
            --limitBAMsortRAM 50000000000 \
            --runThreadN {params.threads} \
            --genomeDir {input.index} \
            --readFilesIn {input.R1} \
            --quantMode TranscriptomeSAM \
            --scoreDelOpen -10000 \
            --scoreInsOpen -10000 \
            --outFilterMultimapNmax 1 \
            --outFilterMismatchNmax 2 \
            --outSAMunmapped Within \
            --readFilesCommand zcat \
            --outFileNamePrefix {params.bam_prefix} \
            --outSAMtype BAM SortedByCoordinate
        '''

##################
# Index Samtools #
##################

rule index_all:
    input:
        expand(str(align_dir / "{sample}.Aligned.sortedByCoord.out.bam.bai"), sample=SAMPLES)


rule samtools_index:
    input:
        aligned = align_dir / "{sample}.Aligned.sortedByCoord.out.bam",
    output:
        index = align_dir / "{sample}.Aligned.sortedByCoord.out.bam.bai",
    shell:
        '''
        samtools index {input.aligned}
        '''


########################
# Quantification QuasR #
########################

rule quantification_input:
    input:
        list(expand(str(align_dir / "{sample}.Aligned.sortedByCoord.out.bam"), sample=SAMPLES)),
    output:
        samplefile = processed_dir / "samplefile.txt"
    run:
        print("samples: " + str(SAMPLES))
        print("filenames: " + str(input))
        df = pd.DataFrame({"FileName": input, "SampleName": SAMPLES})
        df.to_csv(output.samplefile, index=False, sep='\t')

rule quantification:
    input:
        aligned = expand(str(align_dir / "{sample}.Aligned.sortedByCoord.out.bam"), sample=SAMPLES),
        index = expand(str(align_dir / "{sample}.Aligned.sortedByCoord.out.bam.bai"), sample=SAMPLES),
        samplefile = processed_dir / "samplefile.txt"
    params:
        annotfile = str(genome_dir / config["annotation_filename"]),
        genome = str(genome_dir / config["genome_filename"]),
    output:
        count_table = processed_dir / "count_table.csv",
        mapped_unmapped = processed_dir / "mapped_unmapped.csv",
    shell:
        '''
        Rscript {temposeq_quantification_script} {input.samplefile} {params.genome} {params.annotfile} {output.count_table} {output.mapped_unmapped}
        '''




###############
### MultiQC ###
###############


rule multiqc:
    message: "running multiqc for temposeq data"
    input:
        expand(str(align_dir / "{sample}.Aligned.toTranscriptome.out.bam"), sample=SAMPLES),
        expand(str(trim_dir / "{sample}_fastp.json"), sample=SAMPLES),
    output:
        qc_dir / "MultiQC_Report.html"
    log: log_dir / "multiqc.log" 
    shell:
        '''
        multiqc \
        --cl_config "extra_fn_clean_exts: {{ '_fastp.json' }}" \
        --cl_config "sample_names_replace_exact: True" \
        --filename MultiQC_Report.html \
        --interactive \
        --sample-names {metadata_file} \
        -fz {processed_dir}
        mv MultiQC_Report.html MultiQC_Report_data.zip {qc_dir}
        '''
