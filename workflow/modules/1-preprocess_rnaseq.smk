include: "define.smk"
include: "trim.smk"
include: "align.smk"

import hashlib
import json

rsem_genome_filepath = genome_dir / pipeline_config["genome_filename"]
rsem_annotation_filepath = genome_dir / pipeline_config["annotation_filename"]
rsem_index_parameters = {
    "rsem_version": RSEM_version,
    "genome_name": pipeline_config["genome_name"],
}
rsem_index_identity_record = {
    "genome_sha256": sha256_file(rsem_genome_filepath),
    "annotation_sha256": sha256_file(rsem_annotation_filepath),
    "parameters": rsem_index_parameters,
}
rsem_index_identity = hashlib.sha256(
    json.dumps(rsem_index_identity_record, sort_keys=True).encode("utf-8")
).hexdigest()[:16]
RSEM_index = rsem_index_root / f"RSEM_index-{rsem_index_identity}"
RSEM_index_prefix = RSEM_index / pipeline_config["genome_name"]
RSEM_index_manifest = RSEM_index / "index_manifest.txt"

rule pp_rs_all:
    input: 
        processed_dir / "count_table.tsv",
        processed_dir / "isoforms_table.tsv",
        genome_removed

#######################
# QUANTIFICATION RSEM #
#######################

rule RSEM_make_index:
    input:
        genome = rsem_genome_filepath,
        annotations = rsem_annotation_filepath
    params:
        index_dir = RSEM_index,
        index_prefix = RSEM_index_prefix,
        genome_name = pipeline_config["genome_name"],
        manifest = RSEM_index_manifest,
        rsem_version = RSEM_version,
        genome_sha256 = rsem_index_identity_record["genome_sha256"],
        annotation_sha256 = rsem_index_identity_record["annotation_sha256"]
    output:
        RSEM_index_manifest
    conda:
        "../envs/preprocessing.yml"
    benchmark: log_dir / "benchmark.RSEM_make_index.txt"
    shell:
        '''
        actual_rsem_version=$(find "$CONDA_PREFIX/conda-meta" -maxdepth 1 -type f -name 'rsem-*.json' -printf '%f\n' | sed -n 's/^rsem-\([^-]*\)-.*$/\\1/p')
        case "$actual_rsem_version" in
            *"{params.rsem_version}"*) ;;
            *)
                echo "Expected RSEM {params.rsem_version}, but found: $actual_rsem_version" >&2
                exit 1
                ;;
        esac
        mkdir -p {params.index_dir}
        rsem-prepare-reference \
            --gtf {input.annotations} \
            {input.genome} \
            {params.index_prefix}
        printf '%s\n' \
            'genome_sha256={params.genome_sha256}' \
            'annotation_sha256={params.annotation_sha256}' \
            'rsem_version={params.rsem_version}' \
            'genome_name={params.genome_name}' > {params.manifest}
        '''

if pipeline_config["mode"] == "pe":
    rule RSEM:
        input:
            bam = align_dir / "{sample}.Aligned.toTranscriptome.out.bam",
            index = RSEM_index_manifest
        output:
            isoforms = quant_dir / "{sample}.isoforms.results",
            genes = quant_dir / "{sample}.genes.results"
        conda:
            "../envs/preprocessing.yml"
        params:
            threads = pipeline_config["threads"],
            output_prefix =  lambda wildcards : quant_dir / "{}".format(wildcards.sample),
            index_prefix = RSEM_index_prefix
        benchmark: log_dir / "benchmark.{sample}.RSEM_pe.txt"
        threads: pipeline_config["threads"]
        shell:
            '''
            rsem-calculate-expression \
            -p {params.threads} \
            --paired-end \
            --bam {input.bam} \
            --no-bam-output \
            {params.index_prefix} \
            {params.output_prefix}
            '''

if pipeline_config["mode"] == "se":
    rule RSEM:
        input:
            bam = align_dir / "{sample}.Aligned.toTranscriptome.out.bam",
            index = RSEM_index_manifest
        output:
            isoforms = quant_dir / "{sample}.isoforms.results",
            genes = quant_dir / "{sample}.genes.results"
        conda:
            "../envs/preprocessing.yml"
        params:
            threads = pipeline_config["threads"],
            output_prefix =  lambda wildcards : quant_dir / "{}".format(wildcards.sample),
            index_prefix = RSEM_index_prefix
        benchmark: log_dir / "benchmark.{sample}.RSEM_se.txt"
        threads: pipeline_config["threads"]
        shell:
            '''
            rsem-calculate-expression \
            -p {params.threads} \
            --bam {input.bam} \
            --no-bam-output \
            {params.index_prefix} \
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