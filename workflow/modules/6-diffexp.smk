include: "1-define.smk"

rule diff_all:
	input: "reports_complete"

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
        touch("DESeq2_complete")
    conda:
        "../envs/reports.yml"
    benchmark: log_dir / "benchmark.deseq2.txt"
    shell:
        '''
        Rscript scripts/run_DESeq2.R {analysis_folder}
        # rm genome.removed
        '''

######################
### DESeq2 Reports ###
######################

rule deseq_reports:
    message: "generating DESeq2 reports in R..."
    input:
        "DESeq2_complete"
    output:
        touch("reports_complete")
    conda:
        "../envs/reports.yml"
    benchmark: log_dir / "benchmark.deseq_report.txt"
    shell:
        '''
        rm DESeq2_complete
        Rscript scripts/render_DESeq2_report.parallel.R {analysis_folder}
        '''


