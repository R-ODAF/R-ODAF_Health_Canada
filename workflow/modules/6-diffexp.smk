include: "1-define.smk"

rule diff_all:
	input: "analysis_complete"

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
        Rscript scripts/run_DESeq2.R 
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
        Rscript scripts/render_DESeq2_report.parallel.R
        '''

###########################
### Move Deseq2 results ###
###########################

rule move_analysis:
    input: 
        "reports_complete"
    output:
        touch("analysis_complete")
    params:
        analysis_name = deseq_config["analysis_name"]
    shell: 
        '''
        mv output/analysis/most_recent_analysis output/analysis/analysis_{params.analysis_name}_$(date +%Y%m%d_%H%M)
        rm reports_complete
        '''

