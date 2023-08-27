


rule sort:
    input:
        lambda w : 'processing/align/{sample}-{endedness}.sam'.format(
                        sample = w.sample, endedness = 'pe' if config['samples'][w.sample]['is_paired'] else 'se'
                    )    
    output:
        temp('processing/align/{sample}.sorted.bam')
    log:
        'logs/sort/{sample}.sorted.bam.log'
    conda:
        'envs/samtools.yaml'
    resources:
        mem_mb = double_on_failure(config['resources']['sort']['mem_mb']),
        runtime = double_on_failure(config['resources']['sort']['runtime']),
    threads: config['resources']['sort']['threads']
    shell:
        "samtools view -h {input} | samtools sort -O bam -@ {threads} > {output}"



rule markduplicates:
    input:
        rules.sort.output,
    output:
        protected('analysis/samples/{sample}/{sample}.bam'),
    log:
        'logs/markdups/{sample}.bam.log'
    conda:
        'envs/picard.yaml'
    resources:
        mem_mb = double_on_failure(config['resources']['markduplicates']['mem_mb']),
        runtime = double_on_failure(config['resources']['markduplicates']['runtime'])
    threads: config['resources']['markduplicates']['threads']
    shell:
        """
        picard MarkDuplicates \
            I={input} \
            O={output} \
            M={log} && \
        samtools index {output}
        """
