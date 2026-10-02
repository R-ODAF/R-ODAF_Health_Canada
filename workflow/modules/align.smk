import os
import hashlib
import json

# Set STAR parameter based on running on Azure Batch vs locally
if "AZ_BATCH_POOL_ID" in os.environ:
    star_load_mode = "NoSharedMemory"
else:
    star_load_mode = "LoadAndKeep"

################################
### Alignment of reads: STAR ###
################################

if common_config["platform"] =="TempO-Seq":
    STAR_insertion_deletion_penalty = -1000000 # STAR scoring penalty for deletion/insertion, set by biospyder
    STAR_genomeSAindexNbases = 4 # Non-default as specified by BioSpyder
    STAR_multimap_nmax = 1
    STAR_mismatch_nmax = 2
else:
    STAR_insertion_deletion_penalty = -2 # STAR defaults
    STAR_genomeSAindexNbases = 14 # STAR defaults
    STAR_multimap_nmax = 20
    STAR_mismatch_nmax = 999


def sha256_file(path):
    digest = hashlib.sha256()
    with open(path, "rb") as input_file:
        for chunk in iter(lambda: input_file.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


genome_filepath = genome_dir / pipeline_config["genome_filename"]
annotation_filepath = genome_dir / pipeline_config["annotation_filename"]
star_index_parameters = {
    "star_version": STAR_version,
    "sjdbOverhang": 100,
    "genomeSAsparseD": 2,
    "genomeChrBinNbits": 18,
    "genomeSAindexNbases": STAR_genomeSAindexNbases,
    "sjdbGTFfeatureExon": "exon",
}
star_index_identity_record = {
    "genome_sha256": sha256_file(genome_filepath),
    "annotation_sha256": sha256_file(annotation_filepath),
    "parameters": star_index_parameters,
}
star_index_identity = hashlib.sha256(
    json.dumps(star_index_identity_record, sort_keys=True).encode("utf-8")
).hexdigest()[:16]
STAR_index = star_index_root / f"STAR_index-{star_index_identity}"
STAR_index_manifest = STAR_index / "index_manifest.txt"
genome_loaded = sm_temp_dir / f"genome.{star_index_identity}.loaded"
genome_removed = sm_temp_dir / f"genome.{star_index_identity}.removed"

# Build STAR index if not already present
rule STAR_make_index:
    input:
        genome = genome_filepath,
        annotations = annotation_filepath
    params:
        index_dir = STAR_index,
        overhang = star_index_parameters["sjdbOverhang"],
        suffix_array_sparsity = star_index_parameters["genomeSAsparseD"],
        genomeChrBinNbits = star_index_parameters["genomeChrBinNbits"],
        genomeSAindexNbases = star_index_parameters["genomeSAindexNbases"],
        sjdbGTFfeatureExon = star_index_parameters["sjdbGTFfeatureExon"],
        manifest = STAR_index_manifest,
        star_version = star_index_parameters["star_version"],
        genome_sha256 = star_index_identity_record["genome_sha256"],
        annotation_sha256 = star_index_identity_record["annotation_sha256"]
    conda:
        "../envs/preprocessing.yml"
    output:
        STAR_index_manifest
    benchmark: log_dir / "benchmark.STAR_make_index.txt"
    threads: num_threads
    shell:
        '''
        actual_star_version=$(STAR --version 2>&1)
        case "$actual_star_version" in
            *"{params.star_version}"*) ;;
            *)
                echo "Expected STAR {params.star_version}, but found: $actual_star_version" >&2
                exit 1
                ;;
        esac
        STAR \
        --runMode genomeGenerate \
        --genomeDir {params.index_dir} \
        --genomeFastaFiles {input.genome} \
        --sjdbGTFfile {input.annotations} \
        --sjdbOverhang {params.overhang} \
        --runThreadN {threads} \
        --genomeSAsparseD {params.suffix_array_sparsity} \
        --genomeChrBinNbits {params.genomeChrBinNbits} \
        --genomeSAindexNbases {params.genomeSAindexNbases} \
        --sjdbGTFfeatureExon {params.sjdbGTFfeatureExon}
        printf '%s\n' \
            'genome_sha256={params.genome_sha256}' \
            'annotation_sha256={params.annotation_sha256}' \
            'sjdbOverhang={params.overhang}' \
            'genomeSAsparseD={params.suffix_array_sparsity}' \
            'genomeChrBinNbits={params.genomeChrBinNbits}' \
            'genomeSAindexNbases={params.genomeSAindexNbases}' \
            'sjdbGTFfeatureExon={params.sjdbGTFfeatureExon}' \
            'star_version={params.star_version}' > {params.manifest}
        '''

# Check if running on Azure Batch
# If not, load STAR index into shared memory
if "AZ_BATCH_POOL_ID" not in os.environ:
    rule STAR_load:
        input:
            STAR_index_manifest
        output:
            touch(genome_loaded)
        conda:
            "../envs/preprocessing.yml"
        params:
            index = STAR_index
        benchmark: log_dir / "benchmark.STAR_load.txt"
        shell:
            '''
            STAR --genomeLoad LoadAndExit --genomeDir {params.index}
            '''
    # After STAR is run, unload the STAR index from shared memory
    # Delete unnecessary log files made by STAR
    rule STAR_unload:
        input:
            idx = genome_loaded,
            bams = expand(str(align_dir / "{sample}.Aligned.toTranscriptome.out.bam"), sample=SAMPLES)
        output:
            touch(genome_removed)
        conda:
            "../envs/preprocessing.yml"
        params:
            genome_dir = STAR_index
        shell:
            '''
            STAR --genomeLoad Remove --genomeDir {params.genome_dir}
            rm Log.progress.out Log.final.out Log.out SJ.out.tab Aligned.out.sam
            '''


# Run STAR. Depends on config settings.
if pipeline_config["mode"] == "se":
    rule STAR:
        input:
            loaded_index = genome_loaded,
            R1 = trim_dir / "{sample}.fastq.gz"
        output:
            sortedByCoord = align_dir / "{sample}.Aligned.sortedByCoord.out.bam",
            toTranscriptome = align_dir / "{sample}.Aligned.toTranscriptome.out.bam"
        conda:
            "../envs/preprocessing.yml"
        params:
            index = STAR_index,
            penalty = STAR_insertion_deletion_penalty,
            multimap_nmax = STAR_multimap_nmax,
            mismatch_nmax = STAR_mismatch_nmax,
            annotations = genome_dir / pipeline_config["annotation_filename"],
            folder = "{sample}",
            bam_prefix = lambda wildcards : align_dir / "{}.".format(wildcards.sample)
        benchmark: log_dir / "benchmark.{sample}.STAR_se.txt"
        threads: num_threads
        shell:
            '''
            [ -e /tmp/{params.folder} ] && rm -r /tmp/{params.folder}
            STAR \
                --alignEndsType EndToEnd \
                --genomeLoad LoadAndKeep \
                --runThreadN {threads} \
                --genomeDir {params.index} \
                --readFilesIn {input.R1} \
                --quantMode TranscriptomeSAM \
                --limitBAMsortRAM=10737418240 \
                --outTmpDir /tmp/{params.folder} \
                --scoreDelOpen {params.penalty} \
                --scoreInsOpen {params.penalty} \
                --outFilterMultimapNmax {params.multimap_nmax} \
                --outFilterMismatchNmax {params.mismatch_nmax} \
                --readFilesCommand zcat \
                --outFileNamePrefix {params.bam_prefix} \
                --outSAMtype BAM SortedByCoordinate
            '''

if pipeline_config["mode"] == "pe":
    rule STAR:
        input:
            loaded_index = genome_loaded,
            R1 = trim_dir / "{sample}.R1.fastq.gz",
            R2 = trim_dir / "{sample}.R2.fastq.gz"
        output:
            sortedByCoord = align_dir / "{sample}.Aligned.sortedByCoord.out.bam",
            toTranscriptome = align_dir / "{sample}.Aligned.toTranscriptome.out.bam"
        conda:
            "../envs/preprocessing.yml"
        params:
            index = STAR_index,
            annotations = genome_dir / pipeline_config["annotation_filename"],
            folder = "{sample}",
            bam_prefix = lambda wildcards : align_dir / "{}.".format(wildcards.sample)
        resources:
            load=100
        benchmark: log_dir / "benchmark.{sample}.STAR_pe.txt"
        threads: num_threads
        shell:
            '''
            [ -e /tmp/{params.folder} ] && rm -r /tmp/{params.folder}
            STAR \
            --genomeLoad LoadAndKeep \
            --runThreadN {threads} \
            --genomeDir {params.index} \
            --readFilesIn {input.R1} {input.R2} \
            --quantMode TranscriptomeSAM \
            --readFilesCommand zcat \
            --outFileNamePrefix {params.bam_prefix} \
            --outTmpDir /tmp/{params.folder} \
            --outSAMtype BAM SortedByCoordinate
            '''






