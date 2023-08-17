import os
import glob

download_url = '"https://api.ncbi.nlm.nih.gov/datasets/v2alpha/genome/accession/{genome}/download?include_annotation_type=GENOME_FASTA,GENOME_GFF,RNA_FASTA,CDS_FASTA,PROT_FASTA,SEQUENCE_REPORT&filename={genome}.zip"'

rule download_genome:
    output:
        temp('{genome}.zip')
    params:
        url = lambda wildcards : download_url.format(genome = wildcards.genome)
    group: 'download'
    shell:
        'curl -OJX GET {params.url} -H "Accept: application/zip"'

ncbi_fasta_name = 'genomes/{genome}/{genome}_{name}_genomic.fna'


rule move_genome:
    input:
        rules.download_genome.output
    output:
        _dir = directory('genomes/{genome}'),
        gff = 'genomes/{genome}/genomic.gff',
        fasta = 'genomes/{genome}/genomic.fna'
    group: 'download'
    params:
        fasta_download = lambda wildcards : ncbi_fasta_name.format(genome = wildcards.genome, name = "*"),
    shell:
        'mkdir -p genomes/tmp/{wildcards.genome} && mkdir -p genomes/{wildcards.genome} && '
        'unzip {wildcards.genome}.zip -d genomes/tmp/{wildcards.genome} && '\
        'mv genomes/tmp/{wildcards.genome}/ncbi_dataset/data/{wildcards.genome}/* genomes/{wildcards.genome}/ && '
        'mv {params.fasta_download} {output.fasta}'


rule merge_annotations:
    input:
        fastas = expand(rules.move_genome.output.fasta, 
            genome = [metadata['reference'] for species, metadata in config['genomes'].items()]
        ),
        gff = expand(rules.move_genome.output.gff, 
            genome = [metadata['reference'] for species, metadata in config['genomes'].items()]
        )
    output:
        fasta = 'genomes/fasta.fna',
        gff = 'genomes/gff.gff'
    shell:
        "cat {input.fastas} > {output.fasta} && cat {input.gff} | grep -v '#' > {output.gff}"


def get_reference(wildcards):
    if 'reference' in config:
        return config['reference']
    else:
        return rules.merge_annotations.output.fasta


def get_gff(wildcards):
    if 'gff' in config:
        return config['gff']
    else:
        return rules.merge_annotations.output.gff


rule make_index:
    input:
        get_reference
    output:
        'index/bowtie2_index.bt2'
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
        'bowtie2-build {input} {output}'


rule trim_reads_paired:
    input:
        r1 = lambda wildcards: config['samples'][wildcards.sample]['read1'],
        r2 = lambda wildcards: config['samples'][wildcards.sample]['read2']
    output:
        r1 = temp('temp/fastq/{sample}_R1.trim.fastq.gz'),
        r2 = temp('temp/fastq/{sample}_R2.trim.fastq.gz'),
        unmatched1 = temp('temp/{sample}_R1.trim.unpaired.fastq.gz'),
        unmatched2 = temp('temp/{sample}_R2.trim.unpaired.fastq.gz'),
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
        temp('temp/fastq/{sample}.trim.fastq.gz')
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
            MAXINFO:80:0.5 MINLEN:50 AVGQUAL:20 \
            && sleep 10 && du -sh {output} > {log}
        """


rule align_pe:
    input:
        reads = rules.trim_reads_paired.output,
        index = rules.make_index.output
    output:
        temp('temp/align/{sample}-pe.sam')
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
    shell:
        """
        bowtie2 -x {input.index} \
            -1 {input.reads[0]} -2 {input.reads[1]} \
            --threads {threads} \
            --very-sensitive -a --no-unal -S {output} \
            && sleep 10 && du -sh {output} > {log}
        """


rule align_se:
    input:
        reads = rules.trim_reads_unpaired.output,
        index = rules.make_index.output
    output:
        temp('temp/align/{sample}-se.sam')
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
    shell:
        """
        bowtie2 -x {input.index} \
            -U {input.reads} \
            --threads {threads} \
            --very-sensitive -a --no-unal -S {output} \
            && sleep 10 && du -sh {output} > {log}
        """
