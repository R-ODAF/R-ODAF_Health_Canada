include: "define.smk"

# NOTE that preprocessing steps (1-4) must be complete before running this step.
# Smk modules 2-4 are not in "include" statements to allow changing QC params and re-running without Snakemake forcing reruns of previous steps
# Smk module 1 is included because it defines paths used in rules below

# Rule all for if running this module separately from full pipeline
rule qc_all:
	input: qc_dir / "details/samples_removed.txt"


###############
### MultiQC ###
###############

if common_config["platform"] == "DRUG-Seq":
    rule mapstats:
        """
        Extract mapped and unmapped read counts from demultiplexed BAM files
        """
        input: processed_dir / "{library}/demux_bam/{sample}.bam"
        output: qc_dir / "demuxbams_stats/{library}/{sample}.stats.txt"
        conda:
            "../envs/drugseq.yaml"
        shell:
            """
            samtools stats {input} > {output}
            """
    
    rule quantify_dedup:
        """
        Calculate the number of reads in samples when deduplicated vs non-deduplicated.
        Separate output files will be produced for all deduplication methods.
        """
        input:
            dedup = processed_dir / "{library}/{library}_umiDedup-{method}.tsv",
            nodedup = processed_dir / "{library}/{library}_umiDedup-NoDedup.tsv"
        output: qc_dir / "drugseqQC/{library}_umi{method}_dedup_countsums.txt"
        params:
            outdir = qc_dir / "drugseqQC"
        conda:
            "../envs/drugseq.yaml"
        shell:
            """
            Rscript scripts/dedup_stats.R {input.dedup} {input.nodedup} {params.outdir} {wildcards.library} {wildcards.method}
            """
    
    rule fastqc:
        """Quality control of raw FASTQ files"""
        input:
            fastq= os.path.join(raw_dir, "{library}_{read}.fastq.gz")
        output:
            html= qc_dir / "fastqc/{library}_{read}_fastqc.html",
            zip= qc_dir / "fastqc/{library}_{read}_fastqc.zip"
        params:
            outdir= qc_dir / "fastqc"
        conda:
            "../envs/drugseq.yaml"
        resources:
            threads=8,
            mem_mb=1000 
        shell:
            """
            fastqc  --threads {resources.threads} --outdir {params.outdir} {input.fastq}
            """
     
    rule multiqc:
        """
        Consolidate QC files (fastqc, alignment stats) into one report.
        Create summary files used for Studywide QC report
        """
        input:
            files=[f"{qc_dir}/demuxbams_stats/{library}/{sample}.stats.txt" 
                   for library in LIBRARIES 
                   for sample in LIBRARY_SAMPLES[library]],
            fastqc=expand(qc_dir / "fastqc/{library}_{read}_fastqc.zip", 
                        library=LIBRARIES, 
                        read=["R1", "R2"]),
            preprocess_done=sm_temp_dir / "drugseq_preprocess_complete"
        output: 
            qc_dir / "MultiQC_Report.html"
        params:
            outdir= qc_dir,
            demuxstats_dirs = expand(qc_dir / "demuxbams_stats/{library}", library=LIBRARIES)
        conda:
            "../envs/drugseq.yaml"
        shell:
            """
            multiqc \
            --filename MultiQC_Report.html \
            --interactive \
            -fz \
            -o {params.outdir} {params.demuxstats_dirs} output/QC/fastqc/
            """
            
else:
    rule multiqc:
        message: "running multiqc data"
        input:
            expand(str(align_dir / "{sample}.Aligned.toTranscriptome.out.bam"), sample=SAMPLES),
            expand(str(trim_dir / "{sample}_fastp.json"), sample=SAMPLES)
        output:
            qc_dir / "MultiQC_Report.html"
        conda:
            "../envs/preprocessing.yml"
        benchmark: log_dir / "benchmark.multiqc.txt"
        shell:
            '''
            multiqc \
            --cl-config "extra_fn_clean_exts: {{ '_fastp.json' }}" \
            --cl-config "sample_names_replace_exact: True" \
            --filename MultiQC_Report.html \
            --interactive \
            --sample-names {metadata_file} \
            -m fastp -m star -m rsem \
            -fz {raw_dir} {trim_dir} {align_dir} {quant_dir} {processed_dir}
            mv MultiQC_Report.html MultiQC_Report_data.zip {qc_dir}
            '''

#####################
### Study-wide QC ###
#####################

# Set expected inputs for studywide QC, depending on platform
studywide_qc_input = [qc_dir / "MultiQC_Report.html", processed_dir / "count_table.tsv", sm_temp_dir / "genome.removed",]
if common_config["platform"] == "DRUG-Seq":
    # For drugseq, add deduplication count files to input list
    DEDUP_FILES = expand(qc_dir / "drugseqQC/{library}_umi{method}_dedup_countsums.txt", 
                         library=LIBRARIES, method=DEDUP_METHODS_FOR_COMPARISON)
    studywide_qc_input.extend(DEDUP_FILES)


rule studywideqc:
    message: "generating study-wide QC report in R"
    input:
        studywide_qc_input
    output:
        qc_dir / "details/samples_removed.txt"
    conda:
        "../envs/reports.yml"
    benchmark: log_dir / "benchmark.studywide_qc.txt"
    shell:
        '''
        Rscript scripts/render_studywide_QC_report.R
        '''

onerror:
    print("Preprocessing steps must be completed before running this module.")