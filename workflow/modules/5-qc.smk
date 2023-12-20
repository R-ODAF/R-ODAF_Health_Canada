include: "1-define.smk"

# NOTE that preprocessing steps (1-4) must be complete before running this step.
# Smk modules 2-4 are not in "include" statements to allow changing QC params and re-running without Snakemake forcing reruns of previous steps
# Smk module 1 is included because it defines paths used in rules below

# Rule all for if running this module separately from full pipeline
rule qc_all:
	input: qc_dir / "details/samples_removed.txt"


###############
### MultiQC ###
###############

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
        --cl_config "extra_fn_clean_exts: {{ '_fastp.json' }}" \
        --cl_config "sample_names_replace_exact: True" \
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
        "genome.removed"
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