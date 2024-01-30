##################################
### Trimming raw reads : Fastp ###
##################################

if common_config["platform"] =="TempO-Seq":
    length_required_fastp = 50
    trim_tail1_fastp = 1
else:
    length_required_fastp = 36
    trim_tail1_fastp = 0

rule fastp_se:
    input:
        R1 = ancient(str(raw_dir / "{sample}.fastq.gz"))
    output:
        R1 = trim_dir / "{sample}.fastq.gz", #put the pipes back in when done testing!
        json = trim_dir / "{sample}_fastp.json",
        html = trim_dir / "{sample}_fastp.html"
    conda:
        "../envs/preprocessing.yml"
    params:
        front_size = 1,
        front_quality = 3,
        tail_size = 1,
        tail_quality = 3,
        right_size = 4,
        right_quality = 15,
        length_required = length_required_fastp,
        trim_tail1 = trim_tail1_fastp
    benchmark: log_dir / "benchmark.{sample}.fastp_se.txt"
    threads: num_threads
    shell:
        '''
        fastp \
        -i {input.R1} \
        --out1 {output.R1} \
        --json {output.json} \
        --html {output.html} \
        --cut_front \
        --disable_adapter_trimming \
        --cut_front_window_size {params.front_size} \
        --cut_front_mean_quality {params.front_quality} \
        --cut_tail \
        --cut_tail_window_size {params.tail_size} \
        --cut_tail_mean_quality {params.tail_quality} \
        --cut_right \
        --cut_right_window_size {params.right_size} \
        --cut_right_mean_quality {params.right_quality} \
        --length_required {params.length_required} \
        --trim_tail1 {params.trim_tail1} \
        --thread {threads} \
        '''


rule fastp_pe:
    input:
        R1 = ancient(str(raw_dir / "{sample}.R1.fastq.gz")),
        R2 = ancient(str(raw_dir / "{sample}.R2.fastq.gz"))
    output:
        R1 = trim_dir / "{sample}.R1.fastq.gz", #put the pipes back in when done testing!
        R2 = trim_dir / "{sample}.R2.fastq.gz", #put the pipes back in when done testing!
        json = trim_dir / "{sample}_fastp.json",
        html = trim_dir / "{sample}_fastp.html"
    conda:
        "../envs/preprocessing.yml"
    params:
        front_size = 1,
        front_quality = 3,
        tail_size = 1,
        tail_quality = 3,
        right_size = 4,
        right_quality = 15,
        length_required = length_required_fastp
    benchmark: log_dir / "benchmark.{sample}.fastp_pe.txt"
    threads: num_threads
    shell:
        '''
        fastp \
        -i {input.R1} \
        --in2 {input.R2} \
        --out1 {output.R1} \
        --out2 {output.R2} \
        --json {output.json} \
        --html {output.html} \
        --cut_front \
        --cut_front_window_size {params.front_size} \
        --cut_front_mean_quality {params.front_quality} \
        --cut_tail \
        --cut_tail_window_size {params.tail_size} \
        --cut_tail_mean_quality {params.tail_quality} \
        --cut_right \
        --cut_right_window_size {params.right_size} \
        --cut_right_mean_quality {params.right_quality} \
        --length_required {params.length_required} \
        --thread {threads} \
        '''