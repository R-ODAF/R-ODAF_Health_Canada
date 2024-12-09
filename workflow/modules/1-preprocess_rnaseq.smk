include: "define.smk"
include: "trim.smk"
include: "align.smk"

##### Dummy rule all for testing purposes
##### Don't forget to change it!
rule pp_rs_all:
    input: 
        processed_dir / "count_table.tsv",
        processed_dir / "isoforms_table.tsv",
        sm_temp_dir / "genome.removed"

#######################
# QUANTIFICATION RSEM #
#######################

rule RSEM_make_index:
    input:
        genome = genome_dir / pipeline_config["genome_filename"]
    params:
        annotations = genome_dir / pipeline_config["annotation_filename"],
        genome_name = pipeline_config["genome_name"],
        genome_dir = pipeline_config["genomedir"]
    output:
        directory(genome_dir / "RSEM_index")
    conda:
        "../envs/preprocessing.yml"
    benchmark: log_dir / "benchmark.RSEM_make_index.txt"
    shell:
        '''
        mkdir {params.genome_dir}/RSEM_index
        rsem-prepare-reference --gtf {params.annotations} {input.genome} {output}/{params.genome_name}
        '''

if pipeline_config["mode"] == "pe":
    rule RSEM:
        input:
            bam = align_dir / "{sample}.Aligned.toTranscriptome.out.bam",
            index = genome_dir / "RSEM_index"
        output:
            isoforms = quant_dir / "{sample}.isoforms.results",
            genes = quant_dir / "{sample}.genes.results"
        conda:
            "../envs/preprocessing.yml"
        params:
            threads = workflow.cores,
            output_prefix =  lambda wildcards : quant_dir / "{}".format(wildcards.sample),
            index_name = pipeline_config["genome_name"]
        benchmark: log_dir / "benchmark.{sample}.RSEM_pe.txt"
        threads: workflow.cores
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

if pipeline_config["mode"] == "se":
    rule RSEM:
        input:
            bam = align_dir / "{sample}.Aligned.toTranscriptome.out.bam",
            index = genome_dir / "RSEM_index"
        output:
            isoforms = quant_dir / "{sample}.isoforms.results",
            genes = quant_dir / "{sample}.genes.results"
        conda:
            "../envs/preprocessing.yml"
        params:
            threads = workflow.cores,
            output_prefix =  lambda wildcards : quant_dir / "{}".format(wildcards.sample),
            index_name = pipeline_config["genome_name"]
        benchmark: log_dir / "benchmark.{sample}.RSEM_se.txt"
        threads: workflow.cores
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
        genes = processed_dir / "count_table.tsv",
        isoforms = processed_dir / "isoforms_table.tsv"
    conda:
        "../envs/preprocessing.yml"
    benchmark: log_dir / "benchmark.counts_matrix.txt"
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