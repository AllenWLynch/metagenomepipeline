

# This rule merges all the GFF annotations into a single file
# for alignment indexing.
rule merge_annotations:
    input:
        fastas = expand(rules.move_genome.output.fasta, genome = genomes_list),
        gff = expand(rules.move_genome.output.gff, genome = genomes_list),
    output:
        fasta = 'genomes/all/genomic.fa',
        gff = 'genomes/all/genomic.gff',
    shell:
        "cat {input.fastas} > {output.fasta} && " \
        "cat {input.gff} | grep -v '#' | sort -k1,1 -k5,5n > {output.gff}"


# Construct the bowtie index for alignment.
index_prefix = 'processing/index/bowtie2_index' 
rule make_index:
    input:
        get_reference
    output:
        directory('processing/index/')
    resources:
        mem_mb = double_on_failure(config['resources']['index']['mem_mb']),
        runtime = double_on_failure_time(config['resources']['index']['runtime'])
    threads:
        config['resources']['index']['threads']
    conda:
        'envs/bowtie2.yaml'
    log:
        'logs/make_index.log'
    benchmark:
        'benchmark/make_index.tsv'
    params:
        index_prefix = index_prefix
    shell:
        'mkdir -p processing/index && bowtie2-build {input} {params.index_prefix} > {log} 2>&1'


rule test:
    output:
        "test.txt"
    shell:
        "echo HERE > {output}"


rule trim_reads_paired:
    input:
        r1 = lambda wildcards: config['samples'][wildcards.sample]['read1'],
        r2 = lambda wildcards: config['samples'][wildcards.sample]['read2']
    output:
        r1 = temp('processing/fastq/{sample}_R1.trim.fastq.gz'),
        r2 = temp('processing/fastq/{sample}_R2.trim.fastq.gz'),
        unmatched1 = temp('processing/fastq/{sample}_R1.trim.unpaired.fastq.gz'),
        unmatched2 = temp('processing/fastq/{sample}_R2.trim.unpaired.fastq.gz'),
    conda:
        'envs/trimmomatic.yaml'
    resources:
        mem_mb = double_on_failure(config['resources']['trimmomatic']['mem_mb']),
        runtime = double_on_failure_time(config['resources']['trimmomatic']['runtime'])
    threads: config['resources']['trimmomatic']['threads']
    log:
        'logs/trim/{sample}.log'
    benchmark:
        'benchmark/trim/{sample}.tsv'
    params:
        adapters = config['adapters_fa']
    shell:
        """
        trimmomatic \
            PE {input.r1} {input.r2} {output.r1} {output.unmatched1} {output.r2} {output.unmatched2} \
            ILLUMINACLIP:{params.adapters}:2:30:10:8:TRUE \
            MAXINFO:80:0.5 MINLEN:50 AVGQUAL:20 > {log} 2>&1
        """


rule trim_reads_unpaired:
    input:
        lambda wildcards: config['samples'][wildcards.sample]['read1']
    output:
        temp('processing/fastq/{sample}.trim.fastq.gz')
    log:
        'logs/trim/{sample}.log'
    conda:
        'envs/trimmomatic.yaml'
    resources:
        mem_mb = double_on_failure(config['resources']['trimmomatic']['mem_mb']),
        runtime = double_on_failure_time(config['resources']['trimmomatic']['runtime'])
    threads: config['resources']['trimmomatic']['threads']
    message:
        'Trimming reads for SE sample {wildcards.sample}'
    params:
        adapters = config['adapters_fa']
    shell:
        """
        trimmomatic \
            SE {input} {output} \
            ILLUMINACLIP:{params.adapters}:2:30:10 \
            MAXINFO:80:0.5 MINLEN:50 AVGQUAL:20 > {log} 2>&1
        """


rule align_pe:
    input:
        reads = rules.trim_reads_paired.output,
        index = rules.make_index.output
    output:
        sam = 'processing/align/{sample}-pe.sam',
        stats = 'QC/samples/{sample}/flagstat.txt'
    conda:
        'envs/bowtie2.yaml'
    resources:
        mem_mb = double_on_failure(config['resources']['bowtie']['mem_mb']),
        runtime = double_on_failure_time(config['resources']['bowtie']['runtime'])
    threads: config['resources']['bowtie']['threads']
    log:
        'logs/align/{sample}.log'
    benchmark:
        'benchmark/align/{sample}.tsv'
    params:
        index = index_prefix
    shell:
        """
        bowtie2 --rg-id {wildcards.sample} --rg SM:{wildcards.sample} \
            -x {params.index} \
            -1 {input.reads[0]} -2 {input.reads[1]} \
            --threads {threads} \
            --very-sensitive -a --no-unal -S {output.sam} > {log} 2>&1 \
        && samtools flagstat {output.sam} > {output.stats}
        """


rule align_se:
    input:
        reads = rules.trim_reads_unpaired.output,
        index = rules.make_index.output,
    output:
        temp('processing/align/{sample}-se.sam')
    conda:
        'envs/bowtie2.yaml'
    resources:
        mem_mb = double_on_failure(config['resources']['bowtie']['mem_mb']),
        runtime = double_on_failure_time(config['resources']['bowtie']['runtime'])
    threads: config['resources']['bowtie']['threads']
    log:
        'logs/align/{sample}.log'
    benchmark:
        'benchmark/align/{sample}.tsv'
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
