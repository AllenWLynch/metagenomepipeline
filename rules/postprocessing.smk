

rule sort:
    input:
        lambda w : 'processing/align/{sample}-{endedness}.sam'.format(
                        sample = w.sample, endedness = 'pe' if config['samples'][w.sample]['is_paired'] else 'se'
                    )    
    output:
        temp('processing/align/{sample}.sorted.bam')
    conda:
        'envs/samtools.yaml'
    resources:
        mem_mb = double_on_failure(config['resources']['sort']['mem_mb']),
        runtime = double_on_failure_time(config['resources']['sort']['runtime']),
    threads: config['resources']['sort']['threads']
    log:
        'logs/sort/{sample}.log'
    benchmark:
        'benchmark/sort/{sample}.tsv'
    shell:
        "samtools view -h {input} | samtools sort -O bam -@ {threads} > {output} 2> {log}"


rule markduplicates:
    input:
        bams = rules.sort.output,
    output:
        bam = protected('analysis/samples/{sample}/{sample}.bam'),
        metrics = "QC/samples/{sample}/duplicates.txt"
    resources:
        mem_mb = double_on_failure(config['resources']['markduplicates']['mem_mb']),
        runtime = double_on_failure_time(config['resources']['markduplicates']['runtime'])
    threads: config['resources']['markduplicates']['threads']
    log:
        'logs/markdups/{sample}.log'
    benchmark:
        'benchmark/markdups/{sample}.tsv'
    params:
        extra="--CREATE_INDEX true"
    wrapper:
        "v2.6.0/bio/picard/markduplicates"

rule bam_index:
    input:
        rules.markduplicates.output
    output:
        'analysis/samples/{sample}/{sample}.bam.bai'
    conda:
        'envs/samtools.yaml'
    shell:
        "samtools index {input}"


rule get_multimap_stats:
    input:
        bam = rules.markduplicates.output.bam,
        contigs = get_contigs,
    output:
        stats = 'analysis/samples/{sample}/multimap_stats.tsv',
        multimap_sam = 'analysis/samples/{sample}/multimap.sorted-by-name.sam',
    conda:
        "envs/samtools.yaml"
    params:
        scripts = config['_external_scripts']
    shell:
        "bash {params.scripts}/multimap-stats {input.bam} {input.contigs} {output.multimap_sam} {output.stats}"
        
