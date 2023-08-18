
'''
rule filter_sam:
    input:
        lambda w : 'temp/align/{sample}-{endedness}.sam'.format(
                        sample = w.sample, endedness = 'pe' if config['samples'][w.sample]['is_paired'] else 'se'
                    )
    output:
        temp('temp/align/{sample}.filter.sam')
    log:
        'logs/filter/{sample}.filter.sam.log'
    resources:
        mem_mb = double_on_failure(config['resources']['filter']['mem_mb']),
        runtime = double_on_failure(config['resources']['filter']['runtime']),
    threads: 1
    shell:
        """
        python code/filter_sam.py {input} 90 25 > {output} \
            && sleep 10 && du -sh {output} > {log}
        """
'''


rule sort:
    input:
        lambda w : 'temp/align/{sample}-{endedness}.sam'.format(
                        sample = w.sample, endedness = 'pe' if config['samples'][w.sample]['is_paired'] else 'se'
                    )    
    output:
        temp('temp/align/{sample}.sorted.bam')
    log:
        'logs/sort/{sample}.sorted.bam.log'
    conda:
        'envs/samtools.yaml'
    resources:
        mem_mb = double_on_failure(config['resources']['sort']['mem_mb']),
        runtime = double_on_failure(config['resources']['sort']['runtime']),
    threads: config['resources']['sort']['threads']
    shell:
        "samtools view -h {input} | samtools sort -O bam -@ 8 > {output}"

rule markduplicates:
    input:
        rules.sort.output,
    output:
        protected('analysis/samples/{sample}/bam.sorted.mkdup.bam'),
    log:
        'logs/markdups/{sample}.sorted.mkdup.bam.log'
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


rule indexbam:
    input:
        rules.markduplicates.output,
    output:
        protected(rules.markduplicates.output[0] + '.bai')
    conda:
        'envs/samtools.yaml'
    resources:
        mem_mb = double_on_failure(config['resources']['indexbam']['mem_mb']),
        runtime = double_on_failure(config['resources']['indexbam']['runtime'])
    threads: 1
    shell:
        'samtools index {input}'

