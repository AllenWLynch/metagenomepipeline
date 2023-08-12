
rule htseq_count:
    input:
       rules.markduplicates.output,
    output:
        protected('htseq_counts/{sample}.htseq.counts')
    params:
        gtf = config['gtf']
    resources:
        mem_mb = double_on_failure(config['resources']['htseq_count']['mem_mb']),
        runtime = double_on_failure(config['resources']['htseq_count']['runtime'])
    threads: 1
    conda:
        'envs/htseq.yaml'
    shell:
        'htseq-count -r pos -s no -a 0 -t gene {input} {params.gtf} > {output}'


def list_input_samples(w):
    return [rules.markduplicates.output[0].format(sample = sample)
            for sample in config['groups'][w.group]
        ]

rule callvariants:
    input:
        list_input_samples
    output:
        vcf = protected('vcfs/{group}.vcf.gz'),
        stats = protected('vcfs/{group}.stats')
    log:
        'logs/mutect/{group}.mutect.log'
    resources:
        mem_mb = double_on_failure(config['resources']['callvariants']['mem_mb']),
        runtime = double_on_failure(config['resources']['callvariants']['runtime'])
    threads: config['resources']['callvariants']['threads']
    params:
        input_list = lambda w : '-I ' + ' -I '.join(list_input_samples(w)),
        imputed_af = 1/(10*len(config['groups'])),
        reference = config['reference'],
    #conda:
    #    'envs/gatk.yaml'
    shell:
        """
        gatk mutect2 -O {output.vcf} \
            {params.input_list} \
            --af-not-in-resource {params.imputed_af} \
            --reference {params.reference}
        """


rule bamcoverage:
    input:
        bam=rules.markduplicates.output,
        bai=rules.indexbam.output,
    output:
        protected('coverage/{sample}.bamcoverage')
    log:
        "logs/coverage/{sample}.coverage.log",
    resources:
        mem_mb = double_on_failure(config['resources']['bamcoverage']['mem_mb']),
        runtime = double_on_failure(config['resources']['bamcoverage']['runtime'])
    threads: config['resources']['bamcoverage']['threads']
    conda:
        'envs/samtools.yaml'
    shell:
        'samtools depth {input.bam} > {output}'



'''
rule filter_variants:
    input:
        rules.callvariants.output.stats


rule summarize_variants:
    input:
        vcf = rules.callvariants.output.vcf
    output:
        protected('gene_tensors/{sample}.{gene}.variants.pkl')
    shell:
        "echo here"
'''