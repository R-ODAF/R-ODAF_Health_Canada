include: "define.smk"

# STARsolo parameters
cb_start = 1 # Cell barcode start position
cb_len = 14 # Cell barcode length
umi_start = 15 # UMI start position
umi_len = 14 # UMI length
mapq_unique = 60 # Mapping quality threshold
strand = "Forward" # Forward, Reverse, or Unstranded
cell_filter = "None" # STARsolo cell filtering
multimap_nmax = 1 # Maximum number of multimapping alignments to output

rule pp_ds_all:
    input: sm_temp_dir / "drugseq_preprocess_complete"

rule preprocess_output_files:
    input:
        # Demultiplexed BAMs - only valid library-sample combinations
        [f"{processed_dir}/{lib}/demux_bam/{samp}.bam" 
         for lib in LIBRARIES 
         for samp in LIBRARY_SAMPLES[lib]],
        # Per-library count matrices for each dedup method
        expand(processed_dir / "{library}/{library}_umiDedup-{method}.tsv", library=LIBRARIES, method=DEDUP_METHODS),
        # Final combined count table
        processed_dir / "count_table.tsv"
    output:
        preprocess_complete_dummy = sm_temp_dir / "drugseq_preprocess_complete"
    shell:
        """
        touch {output.preprocess_complete_dummy}
        """

rule whitelist_barcodes:
    """Format all barcodes for both STARsolo and Picard demultiplexing
    
    This single rule reads the metadata file
    and outputs a single barcode whitelist for STARsolo (one barcode per line, NO HEADER)
    
    This is a project-wide file used for all libraries.
    """
    input:
        metadata=metadata_file
    output:
        starsolo= output_dir / "barcodes/barcode_whitelist.txt",
    params:
        sample_id_col=config["pipeline"]["sample_id"],
        barcode_col=config["pipeline"]["sample_barcode"]
    run:
        import pandas as pd
        
        # Load metadata
        df = pd.read_csv(input.metadata, sep="\t")
        
        # Output STARsolo whitelist (one barcode per line, unique barcodes, NO HEADER)
        barcodes = df[params.barcode_col].unique().tolist()
        
        with open(output.starsolo, 'w') as f:
            for barcode in barcodes:
                f.write(f"{barcode}\n")
        
        print(f"✓ Created STARsolo barcode file: {output.starsolo} ({len(barcodes)} unique barcodes)")


rule barcodes_w_IDs:
    """Format barcodes for identifying STARsolo output and demultiplexing with Picard 
    
    This rule outputs a tsv file with barcodes and samleIDs for each library.
    """
    input:
        metadata=metadata_file
    output:
        picard= output_dir / "barcodes/{library}_barcodes_sampleIDs.txt"
    params:
        sample_id_col=config["pipeline"]["sample_id"],
        barcode_col=config["pipeline"]["sample_barcode"]
    run:
        import pandas as pd
        
        # Load metadata
        df = pd.read_csv(input.metadata, sep="\t")
        
        # Output demux info (WITH HEADER: sample_id and barcode)
        # libraries = df['library_ID'].unique().tolist()
        
        dffilt = df[df['library_ID'].str.contains(wildcards.library)]

        picard_data = []
        for _, row in dffilt.iterrows():
            sample_id = row[params.sample_id_col]
            barcode = row[params.barcode_col]
            picard_data.append({'sample_id': sample_id, 'barcode': barcode})
        
        # Remove duplicates (in case same sample_id/barcode used across libraries)
        picard_df = pd.DataFrame(picard_data).drop_duplicates()
         
        # Write with header
        picard_df.to_csv(output.picard, sep='\t', index=False, header=True)
      
        print(f"✓ Created Picard demux info: {output.picard} ({len(picard_df)} unique sample_id/barcode combinations)")


####################################################################################
# STARsolo alignment and matrix production
####################################################################################

if pipeline_config["include_ercc"] == True:
    rule combine_reference_genomes:
        """Combine reference genome and annotation with ERCC spike-ins if requested
        
        Output is stored in genomedir so it can be reused across projects.
        Only runs if include_ercc is TRUE (triggered by rule star_index)
        """
        input:
            primary_fasta=genome_filepath,
            primary_gtf=annotation_filepath,
            ercc_fasta=ercc_fasta,
            ercc_gtf=ercc_gtf
        output:
            combined_fasta=genome_fasta,
            combined_gtf=genome_gtf
        shell:
            """
            cat {input.primary_fasta} {input.ercc_fasta} > {output.combined_fasta}
            cat {input.primary_gtf} {input.ercc_gtf} > {output.combined_gtf}
            """


rule star_index:
    """Create STAR genome index including ERCC spike-ins if applicable
    
    Index is stored in genomedir so it can be reused across projects.
    Uses combined files if include_ercc is TRUE, otherwise uses primary files.
    """
    input:
        fasta=genome_fasta,
        gtf=genome_gtf
    output:
        genome=index_genome,
    params:
        genomeDir=index_dir,
        threads=workflow.cores
    conda:
        "../envs/drugseq.yaml"
    resources:
        mem_mb=40000,
        threads=workflow.cores
    shell:
        """
        STAR --runMode genomeGenerate \
            --genomeDir {params.genomeDir} \
            --genomeFastaFiles {input.fasta} \
            --sjdbGTFfile {input.gtf} \
            --runThreadN {params.threads}
        """

# rule STAR_load:
#         input:
#             genome_dir / "STAR_index"
#         output:
#             touch(sm_temp_dir / "genome.loaded")
#         conda:
#             "../envs/preprocessing.yml"
#         params:
#             index = index_dir
#         benchmark: log_dir / "benchmark.STAR_load.txt"
#         shell:
#             '''
#             STAR --genomeLoad LoadAndExit --genomeDir {params.index}
#             '''

rule STAR_unload:
    input:
        # idx = sm_temp_dir / "genome.loaded",
        bams = expand(processed_dir / "{library}/STARsolo/Aligned.sortedByCoord.out.bam", library=LIBRARIES)
    output:
        touch(sm_temp_dir / "genome.removed")
    conda:
        "../envs/drugseq.yaml"
    params:
        genome_dir = index_dir
    shell:
        """
        STAR --genomeLoad Remove --genomeDir {params.genome_dir}
        """

rule starsolo:
    """Align reads and generate count matrices using STARsolo"""
    input:
        r2=os.path.join(raw_dir, "{library}_R2.fastq.gz"),
        r1=os.path.join(raw_dir, "{library}_R1.fastq.gz"),
        genome=index_genome,
        barcodes=output_dir / "barcodes/barcode_whitelist.txt"
    output:
        processed_dir / "{library}/STARsolo/Aligned.sortedByCoord.out.bam",
        processed_dir / "{library}/STARsolo/Solo.out/Gene/raw/features.tsv",
        processed_dir / "{library}/STARsolo/Solo.out/Gene/raw/barcodes.tsv",
        processed_dir / "{library}/STARsolo/Solo.out/Barcodes.stats",
        processed_dir / "{library}/STARsolo/Log.final.out",
        # MTX files for each dedup method
        *[f"{processed_dir}/{{library}}/STARsolo/Solo.out/Gene/raw/umiDedup-{method}.mtx" 
          for method in DEDUP_METHODS]
    params:
        outprefix= "output/processed/{library}/STARsolo/", # this path MUST be hardcoded for {library} wildcard to work. I don't understand why
        genomeDir=index_dir,
        threads=workflow.cores,
        cb_start=cb_start,
        cb_len=cb_len,
        umi_start=umi_start,
        umi_len=umi_len,
        mapq=mapq_unique,
        bam_threads=pipeline_config["threads"],
        strand=strand,
        umi_dedup=pipeline_config["umi_dedup_method"],
        cell_filter=cell_filter,
        multimap_nmax=multimap_nmax,
        max_RAM=pipeline_config["max_memory"]
    conda:
        "../envs/drugseq.yaml"
    threads: workflow.cores # Run a single starsolo job at a time
    shell:
        """
        STAR --runMode alignReads \
            --genomeLoad LoadAndKeep \
            --limitBAMsortRAM={params.max_RAM} \
            --outSAMmapqUnique {params.mapq} \
            --runThreadN {params.threads} \
            --outSAMunmapped Within \
            --soloStrand {params.strand} \
            --quantMode GeneCounts \
            --outBAMsortingThreadN {params.bam_threads} \
            --genomeDir {params.genomeDir} \
            --soloType CB_UMI_Simple \
            --soloCBstart {params.cb_start} \
            --soloCBlen {params.cb_len} \
            --soloUMIstart {params.umi_start} \
            --soloUMIlen {params.umi_len} \
            --soloUMIdedup {params.umi_dedup} \
            --soloCellFilter {params.cell_filter} \
            --soloCBwhitelist {input.barcodes} \
            --soloBarcodeReadLength 0 \
            --soloFeatures Gene \
            --outSAMattributes NH HI nM AS CR UR CB UB GX GN sS sQ sM \
            --outFilterMultimapNmax {params.multimap_nmax} \
            --readFilesCommand zcat \
            --outSAMtype BAM SortedByCoordinate \
            --outFileNamePrefix {params.outprefix} \
            --readFilesIn {input.r2} {input.r1}
        """

rule mtx_to_counts:
    """Convert STARsolo MTX format to tsv count matrix"""
    input:
        mtx= processed_dir / "{library}/STARsolo/Solo.out/Gene/raw/umiDedup-{method}.mtx",
        # features=processed_dir / "{library}/STARsolo/Solo.out/Gene/raw/features.tsv",
        barcodes= output_dir / "barcodes/{library}_barcodes_sampleIDs.txt"
    output:
        counts= processed_dir / "{library}/{library}_umiDedup-{method}.tsv"
    params:
        matrix_dir= processed_dir / "{library}/STARsolo/Solo.out/Gene/raw"
    conda:
        "../envs/drugseq.yaml"
    shell:
        "Rscript scripts/mtx_to_counts.R {input.mtx} {wildcards.library} {input.barcodes} {output.counts}"


rule combine_counttables:
    """
    Combine count tables from multiple libraries into a single project-wide count table for downstream analysis.
    Combines tables deduplicated with the definitive UMI deduplication method specified in the config file.
    """
    input: 
        tables=expand(processed_dir / f"{{library}}/{{library}}_umiDedup-{dedup_method_for_downstream}.tsv", library = LIBRARIES),
        metadata=metadata_file
    output: 
        counttable = processed_dir / "count_table.tsv", # To fit into current R-ODAF needs
        # dummy= sm_temp_dir / "genome.removed" # Super annoying, but other preprocessing paths make it, needed for QC to run
    shell:
        """
        python scripts/combine_counts.py {input.metadata} {input.tables} {output.counttable}
        """


####################################################################################
# Demultiplex BAM files and get mapping stats
####################################################################################

rule picard_demultiplex_sample:
    """
    Demultiplex a single BAM file by sample barcode using Picard.
    Uses CB (corrected barcode), not CR (raw barcode) as in Alithea manual
    """
    input:
        bam=processed_dir / "{library}/STARsolo/Aligned.sortedByCoord.out.bam",
        demux_info=output_dir / "barcodes/{library}_barcodes_sampleIDs.txt"
    output:
        bam=processed_dir / "{library}/demux_bam/{sample}.bam",
    params:
        tag_value=lambda wildcards: get_barcode_for_sample(wildcards)
    conda:
        "../envs/drugseq.yaml"
    shell:
        """
        picard FilterSamReads \
            I={input.bam} \
            O={output.bam} \
            TAG=CB TAG_VALUE={params.tag_value} \
            FILTER=includeTagValues
        """

