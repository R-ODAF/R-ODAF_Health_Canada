include: "define.smk"

# NOTE that preprocessing steps (1-4) and QC steps (5) must be complete before running this module.
# Smk modules 2-5 are not in "include" statements to allow changing deseq2 params and re-running without Snakemake forcing reruns of previous steps
# Smk module 1 is included because it defines paths used in rules below

# Rule all for if running this module separately from full pipeline
rule diff_all:
	input: dummy_dir / "reports_complete"

# Set up analysis folder name with timestamp
from datetime import datetime
t = datetime.now()

analysis_folder = "analysis_" + deseq_config["analysis_name"] + "_" + t.strftime('%Y%m%d-%H%M')

##############
### DESeq2 ###
##############

rule deseq2:
    message: "running DESeq2..."
    input:
        qc_dir / "details/samples_removed.txt"
    output:
        touch(dummy_dir / "DESeq2_complete")
    conda:
        "../envs/reports.yml"
    benchmark: log_dir / "benchmark.deseq2.txt"
    shell:
        '''
        Rscript scripts/run_DESeq2.R {analysis_folder}
        '''

######################
### DESeq2 Reports ###
######################

rule deseq_reports:
    message: "generating DESeq2 reports in R..."
    input:
        dummy_dir / "DESeq2_complete"
    output:
        touch(dummy_dir / "reports_complete")
    conda:
        "../envs/reports.yml"
    benchmark: log_dir / "benchmark.deseq_report.txt"
    shell:
        '''
        Rscript scripts/render_DESeq2_report.parallel.R {analysis_folder}

        conda env export > output/analysis/{analysis_folder}/Pipeline_record/reports_env.yml
        '''

onerror:
    print("Preprocessing and QC steps must be completed before running this module.")
