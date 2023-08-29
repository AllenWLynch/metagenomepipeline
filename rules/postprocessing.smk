


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
        rules.sort.output,
    output:
        protected('analysis/samples/{sample}/{sample}.bam'),
    conda:
        'envs/picard.yaml'
    resources:
        mem_mb = double_on_failure(config['resources']['markduplicates']['mem_mb']),
        runtime = double_on_failure_time(config['resources']['markduplicates']['runtime'])
    threads: config['resources']['markduplicates']['threads']
    log:
        'logs/markdups/{sample}.log'
    benchmark:
        'benchmark/markdups/{sample}.tsv'
    shell:
        """
        picard MarkDuplicates \
            I={input} \
            O={output} \
            M={log} > {log} 2>&1 && \
        samtools index {output}
        """

rule index:
    input:
        rules.markduplicates.output
    output:
        'analysis/samples/{sample}/{sample}.bam.bai'
    conda:
        'envs/samtools.yaml'
    shell:
        "samtools index {input}"
