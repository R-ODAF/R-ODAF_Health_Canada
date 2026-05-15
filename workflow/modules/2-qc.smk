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
        input: "output/{library}/demux_bam/{sample}.bam"
        output: "output/QC/demuxbams_stats/{library}/{sample}.stats.txt"
        conda:
            "../envs/drugseq.yaml"
        shell:
            """
            samtools stats {input} > {output}
            """
     
    rule multiqc:
        """
        Consolidate QC files (fastqc, alignment stats) into one report.
        Create summary files used for Studywide QC report
        """
        input:
            files=expand("output/QC/demuxbams_stats/{library}/{sample}.stats.txt", 
                        library=LIBRARIES, 
                        sample=[s for lib in LIBRARIES for s in LIBRARY_SAMPLES[lib]]),
            fastqc=expand("output/QC/fastqc/{library}_{read}_fastqc.zip", 
                        library=LIBRARIES, 
                        read=["R1", "R2"])
        output: 
            qc_dir / "MultiQC_Report.html"
        params:
            outdir="output/QC",
            # This works for projects with a single library, I suspect it will break for multiple libraries :(
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

rule studywideqc:
    message: "generating study-wide QC report in R"
    input:
        qc_dir / "MultiQC_Report.html",
        processed_dir / "count_table.tsv",
        sm_temp_dir / "genome.removed"
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