

# Construct the bowtie index for alignment.
# 
index_prefix = 'index/bowtie2_index.bt2'
rule make_index:
    input:
        get_reference
    output:
        index_prefix + '.1.bt2'
    log:
        'logs/index/index.log'
    resources:
        mem_mb = double_on_failure(config['resources']['index']['mem_mb']),
        runtime = double_on_failure(config['resources']['index']['runtime'])
    threads:
        config['resources']['index']['threads']
    message:
        'Constructing bowtie index'
    conda:
        'envs/bowtie2.yaml'
    shell:
        'bowtie2-build {input} ' + index_prefix


rule trim_reads_paired:
    input:
        r1 = lambda wildcards: config['samples'][wildcards.sample]['read1'],
        r2 = lambda wildcards: config['samples'][wildcards.sample]['read2']
    output:
        r1 = temp('processing/fastq/{sample}_R1.trim.fastq.gz'),
        r2 = temp('processing/fastq/{sample}_R2.trim.fastq.gz'),
        unmatched1 = temp('processing/fastq/{sample}_R1.trim.unpaired.fastq.gz'),
        unmatched2 = temp('processing/fastq/{sample}_R2.trim.unpaired.fastq.gz'),
    log:
        'logs/trim/{sample}.trim.log'
    conda:
        'envs/trimmomatic.yaml'
    resources:
        mem_mb = double_on_failure(config['resources']['trimmomatic']['mem_mb']),
        runtime = double_on_failure(config['resources']['trimmomatic']['runtime'])
    threads: config['resources']['trimmomatic']['threads']
    message:
        'Trimming reads for PE sample {wildcards.sample}'
    params:
        adapters = config['adapter_fa']
    shell:
        """
        trimmomatic \
            PE {input.r1} {input.r2} {output.r1} {output.unmatched1} {output.r2} {output.unmatched2} \
            ILLUMINACLIP:{params.adapters}:2:30:10:8:TRUE \
            MAXINFO:80:0.5 MINLEN:50 AVGQUAL:20
        """


rule trim_reads_unpaired:
    input:
        lambda wildcards: config['samples'][wildcards.sample]['read1']
    output:
        temp('processing/fastq/{sample}.trim.fastq.gz')
    log:
        'logs/trim/{sample}.trim.log'
    conda:
        'envs/trimmomatic.yaml'
    resources:
        mem_mb = double_on_failure(config['resources']['trimmomatic']['mem_mb']),
        runtime = double_on_failure(config['resources']['trimmomatic']['runtime'])
    threads: config['resources']['trimmomatic']['threads']
    message:
        'Trimming reads for SE sample {wildcards.sample}'
    params:
        adapters = config['adapter_fa']
    shell:
        """
        trimmomatic \
            SE {input} {output} \
            ILLUMINACLIP:{params.adapters}:2:30:10 \
            MAXINFO:80:0.5 MINLEN:50 AVGQUAL:20
        """


rule align_pe:
    input:
        reads = rules.trim_reads_paired.output,
        index = rules.make_index.output
    output:
        sam = 'processing/align/{sample}-pe.sam',
        stats = 'QC/samples/{sample}/flagstat.txt'
    log:
        'logs/align/{sample}.sam.log'
    conda:
        'envs/bowtie2.yaml'
    resources:
        mem_mb = double_on_failure(config['resources']['bowtie']['mem_mb']),
        runtime = double_on_failure(config['resources']['bowtie']['runtime'])
    threads: config['resources']['bowtie']['threads']
    message:
        'Aligning reads for {wildcards.sample}'
    params:
        index = index_prefix
    shell:
        """
        bowtie2 --rg-id {wildcards.sample} --rg SM:{wildcards.sample} \
            -x {params.index} \
            -1 {input.reads[0]} -2 {input.reads[1]} \
            --threads {threads} \
            --very-sensitive -a --no-unal -S {output.sam} 2> {log} \
        && samtools flagstat {output.sam} > {output.stats}
        """


rule align_se:
    input:
        reads = rules.trim_reads_unpaired.output,
        index = rules.make_index.output,
    output:
        temp('processing/align/{sample}-se.sam')
    log:
        'logs/align/{sample}.sam.log'
    conda:
        'envs/bowtie2.yaml'
    resources:
        mem_mb = double_on_failure(config['resources']['bowtie']['mem_mb']),
        runtime = double_on_failure(config['resources']['bowtie']['runtime'])
    threads: config['resources']['bowtie']['threads']
    message:
        'Aligning reads for {wildcards.sample}'
    params:
        index = index_prefix
    shell:
        """
        bowtie2 -x {params.index} \
            -U {input.reads} \
            --threads {threads} \
            --very-sensitive -a --no-unal -S {output} \
            && sleep 10 && du -sh {output} > {log}
        """
