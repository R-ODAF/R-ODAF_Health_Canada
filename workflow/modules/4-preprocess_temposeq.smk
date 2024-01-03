include: "1-define.smk"
include: "2-trim.smk"
include: "3-align.smk"


# Rule all for if only running pipeline to this step
rule pp_ts_all:
    input: 
        processed_dir / "count_table.tsv",
        "genome.removed"

########################
# Quantification QuasR #
########################

rule samtools_index:
    input:
        aligned = align_dir / "{sample}.Aligned.sortedByCoord.out.bam"
    output:
        index = align_dir / "{sample}.Aligned.sortedByCoord.out.bam.bai"
    conda:
        "../envs/temposeqr.yml"
    benchmark: log_dir / "benchmark.{sample}.samtools_index.txt"
    shell:
        '''
        samtools index {input.aligned}
        '''

rule quantification_input:
    input:
        list(expand(str(align_dir / "{sample}.Aligned.sortedByCoord.out.bam"), sample=SAMPLES))
    output:
        samplefile = processed_dir / "samplefile.txt"
    benchmark: log_dir / "benchmark.quantification_input.txt"
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
        annotfile = str(genome_dir / pipeline_config["annotation_filename"]),
        threads = workflow.cores,
        genome = str(genome_dir / pipeline_config["genome_filename"]),
        temposeq_quantification_script = str(main_dir / "scripts/temposeq/temposeq_quantification.R")
    threads: workflow.cores
    output:
        count_table = processed_dir / "count_table.tsv"
    conda:
        "../envs/temposeqr.yml"
    benchmark: log_dir / "benchmark.quantification.txt"
    shell:
        '''
        Rscript \
        {params.temposeq_quantification_script} \
        {input.samplefile} \
        {params.genome} \
        {params.annotfile} \
        {output.count_table} \
        {params.threads}
        '''