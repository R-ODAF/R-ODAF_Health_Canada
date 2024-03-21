import os

# Set STAR parameter based on running on Azure Batch vs locally
if "AZ_BATCH_POOL_ID" in os.environ:
    star_load_mode = "NoSharedMemory"
else:
    star_load_mode = "LoadAndKeep"

################################
### Alignment of reads: STAR ###
################################

if common_config["platform"] =="TempO-Seq":
    STAR_index = genome_dir / "STAR_index"
    STAR_insertion_deletion_penalty = -1000000 # STAR scoring penalty for deletion/insertion, set by biospyder
    STAR_genomeSAindexNbases = 4 # Non-default as specified by BioSpyder
    STAR_multimap_nmax = 1
    STAR_mismatch_nmax = 2
else:
    STAR_index = genome_dir / "STAR_index"
    STAR_insertion_deletion_penalty = -2 # STAR defaults
    STAR_genomeSAindexNbases = 14 # STAR defaults
    STAR_multimap_nmax = 20
    STAR_mismatch_nmax = 999

# Build STAR index if not already present
rule STAR_make_index:
    input:
        genome = ancient(genome_dir / pipeline_config["genome_filename"])
    params:
        index_dir = STAR_index,
        annotations = genome_dir / pipeline_config["annotation_filename"],
        overhang = 100,
        suffix_array_sparsity = 2, # bigger for smaller (RAM), slower indexing
        genomeChrBinNbits = 18, # might need to mess with this for different genomes - default is 18
        genomeSAindexNbases = STAR_genomeSAindexNbases,
        sjdbGTFfeatureExon = "exon" # STAR default
    conda:
        "../envs/preprocessing.yml"
    output:
        directory(genome_dir / "STAR_index")
    benchmark: log_dir / "benchmark.STAR_make_index.txt"
    threads: num_threads
    shell:
        '''
        STAR \
        --runMode genomeGenerate \
        --genomeDir {params.index_dir} \
        --genomeFastaFiles {input.genome} \
        --sjdbGTFfile {params.annotations} \
        --sjdbOverhang {params.overhang} \
        --runThreadN {threads} \
        --genomeSAsparseD {params.suffix_array_sparsity} \
        --genomeChrBinNbits {params.genomeChrBinNbits} \
        --genomeSAindexNbases {params.genomeSAindexNbases} \
        --sjdbGTFfeatureExon {params.sjdbGTFfeatureExon}
        '''

# Check if running on Azure Batch
# If not, load STAR index into shared memory
if "AZ_BATCH_POOL_ID" not in os.environ:
    rule STAR_load:
        input:
            genome_dir / "STAR_index"
        output:
            touch("genome.loaded")
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
            idx = "genome.loaded",
            bams = expand(str(align_dir / "{sample}.Aligned.toTranscriptome.out.bam"), sample=SAMPLES)
        output:
            touch("genome.removed")
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
            loaded_index = "genome.loaded",
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
            alignment_dir = align_dir,
            bam_prefix = lambda wildcards : align_dir / "{}.".format(wildcards.sample),
            load_mode = star_load_mode
        benchmark: log_dir / "benchmark.{sample}.STAR_pe.txt"
        threads: num_threads
        shell:
            '''
            if test -d "{params.alignment_dir}"; then
                echo "Directory {params.alignment_dir} exists."
            else
                echo "Directory {params.alignment_dir} does not exist."
            fi

            [ -e /tmp/{params.folder} ] && rm -r /tmp/{params.folder}
            STAR \
                --alignEndsType EndToEnd \
                --genomeLoad {params.load_mode} \
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
            loaded_index = "genome.loaded",
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
            bam_prefix = lambda wildcards : align_dir / "{}.".format(wildcards.sample),
            load_mode = star_load_mode
        resources:
            load=100
        benchmark: log_dir / "benchmark.{sample}.STAR_se.txt"
        threads: num_threads
        shell:
            '''
            [ -e /tmp/{params.folder} ] && rm -r /tmp/{params.folder}

            

            STAR \
            --genomeLoad {params.load_mode} \
            --runThreadN {threads} \
            --genomeDir {params.index} \
            --readFilesIn {input.R1} {input.R2} \
            --quantMode TranscriptomeSAM \
            --readFilesCommand zcat \
            --outFileNamePrefix {params.bam_prefix} \
            --outTmpDir /tmp/{params.folder} \
            --outSAMtype BAM SortedByCoordinate
            '''






