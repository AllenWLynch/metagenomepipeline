

# java -jar /broad/smillie-data/home/akumbhar/bin/Trimmomatic-0.39/trimmomatic-0.39.jar

def double_on_failure(base_resources):
    def _double_on_failure(wildcards, attempt):
        return base_resources * attempt

    return _double_on_failure


rule make_index:
    input:
        config['reference']
    output:
        'bowtie2_index.bt2'
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
    shell:
        """
        trimmomatic \
            PE {input.r1} {input.r2} {output.r1} {output.unmatched1} {output.r2} {output.unmatched2} \
            ILLUMINACLIP:/broad/smillie-data/db/adapters.mgx.fa:2:30:10:8:TRUE \
            MAXINFO:80:0.5 MINLEN:50 AVGQUAL:20 \
            && sleep 10 && du -sh {output.r1} > {log}
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
    shell:
        """
        trimmomatic \
            SE {input} {output} \
            ILLUMINACLIP:/broad/smillie-data/db/adapters.mgx.fa:2:30:10 \
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
