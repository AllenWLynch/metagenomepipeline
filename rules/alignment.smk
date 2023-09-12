

# This rule merges all the GFF annotations into a single file
# for alignment indexing.
rule merge_annotations:
    input:
        fastas = expand(rules.unpack_genome.output.fasta, genome = genomes_list),
        gff = expand(rules.unpack_genome.output.gff, genome = genomes_list),
        contigs = expand(rules.make_contigs_file.output, genome = genomes_list),
    output:
        fasta = 'genomes/all/genomic.fa',
        gff = 'genomes/all/genomic.gff',
        index = 'genomes/all/genomic.fa.fai',
        chromsizes = 'genomes/all/genomic.chromsizes',
        dict = 'genomes/all/genomic.dict',
        contigs = 'genomes/all/contigs.tsv'
    conda:
        "envs/samtools.yaml"
    shell:
        """
        cat {input.fastas} > {output.fasta} && \
        cat {input.gff} | grep -v '#' | sort -k1,1 -k5,5n > {output.gff} && \
        \
        samtools faidx {output.fasta} && \
        cut -f1,2 {output.fasta}.fai > {output.chromsizes} && \
        samtools dict {output.fasta} -o {output.dict} && \
        \
        head -n1 {input.contigs[0]} > {output.contigs}
        cat {input.contigs} | grep -v '^#genome' >> {output.contigs}
        """


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



def get_reads_paired(wildcards):
    if 'is_trimmed' in config['samples'][wildcards.sample] and \
        config['samples'][wildcards.sample]['is_trimmed']:

        return [
            config['samples'][wildcards.sample]['read1'],
            config['samples'][wildcards.sample]['read2']
        ]
    else:
        return [
            rules.trim_reads_paired.output.r1,
            rules.trim_reads_paired.output.r2
        ]


rule align_pe:
    input:
        reads = get_reads_paired,
        index = rules.make_index.output
    output:
        sam = temp('processing/align/{sample}-pe.sam'),
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


def get_reads_se(wildcards):
    
    if 'is_trimmed' in config['samples'][wildcards.sample] and \
        config['samples'][wildcards.sample]['is_trimmed']:

        return config['samples'][wildcards.sample]['read1'],

    else:
        return rules.trim_reads_unpaired.output.r1


rule align_se:
    input:
        reads = rules.trim_reads_unpaired.output,
        index = rules.make_index.output,
    output:
        sam = temp('processing/align/{sample}-se.sam'),
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
        bowtie2 -x {params.index} \
            -U {input.reads} \
            --threads {threads} \
            --very-sensitive -a --no-unal -S {output} > {log} 2>&1 \
        && samtools flagstat {output.sam} > {output.stats}
        """
