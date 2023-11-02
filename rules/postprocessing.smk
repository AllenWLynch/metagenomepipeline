

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
        "samtools view -h {input} | samtools sort -O bam -@ {threads} -T $TMPDIR > {output} 2> {log}"


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
    wrapper:
        "v2.6.0/bio/picard/markduplicates"


rule bam_index:
    input:
        rules.markduplicates.output.bam
    output:
        'analysis/samples/{sample}/{sample}.bam.bai'
    conda:
        'envs/samtools.yaml'
    shell:
        "samtools index {input}"


rule bamcoverage:
    input:
        bam = get_bam,
        index = get_bai,
        chromsizes = get_chromsizes,
    output:
        bedgraph = temp('analysis/samples/{sample}/coverage.bedgraph'),
        bigwig = 'analysis/samples/{sample}/coverage.bigwig'
    resources:
        mem_mb = double_on_failure(config['resources']['bamcoverage']['mem_mb']),
        runtime = double_on_failure_time(config['resources']['bamcoverage']['runtime'])
    threads: config['resources']['bamcoverage']['threads']
    conda:
        'envs/bedtools.yaml'
    log:
        'logs/coverage/{sample}.log',
    benchmark:
        'benchmark/coverage/{sample}.tsv'
    params:
        filter_str = '-b -q {quality} -F 0x400 -F 0x100 -F 0x800'.format(quality=config['min_count_quality'])
    shell:
        """
        samtools view {params.filter_str} {input.bam} -h | \
        bedtools genomecov -ibam - -bga | sort -k1,1 -k2,2n > {output.bedgraph} 2> {log} && \
        bedGraphToBigWig {output.bedgraph} {input.chromsizes} {output.bigwig}
        """


use rule bamcoverage as pileup_primaries with:
    output:
        bedgraph = temp('analysis/samples/{sample}/multimapping/primary-coverage.bedgraph'),
        bigwig = 'analysis/samples/{sample}/multimapping/primary-coverage.bigwig',
    params:
        filter_str = '-b -q 0 -F -0x800 -F 0x400 -F 0x100'


rule get_multimap_stats:
    input:
        bam = rules.markduplicates.output.bam,
        contigs = get_contigs,
    output:
        multimap_sam = 'analysis/samples/{sample}/multimapping/multimap.sorted-by-name.sam',
        multimap_bam = 'analysis/samples/{sample}/multimapping/multimap.bam'
    conda:
        "envs/samtools.yaml"
    params:
        scripts = config['_external_scripts']
    resources:
        mem_mb = double_on_failure(config['resources']['sort']['mem_mb']),
        runtime = double_on_failure_time(config['resources']['sort']['runtime']),
    threads: config['resources']['sort']['threads']
    shell:
        "bash {params.scripts}/multimap-stats {input.bam} {input.contigs} {output.multimap_sam} && "
        "samtools sort -O bam -@ {threads} -T $TMPDIR {output.multimap_sam} > {output.multimap_bam}"


use rule bamcoverage as pileup_multimapping with:
    input:
        bam = rules.get_multimap_stats.output.multimap_bam,
        chromsizes = get_chromsizes,
    output:
        bedgraph = temp('analysis/samples/{sample}/multimapping/multimap-coverage.bedgraph'),
        bigwig = 'analysis/samples/{sample}/multimapping/multimap-coverage.bigwig',
    params:
        filter_str = '-b -F 0x100'


rule aggregate_multimapping:
    input:
        bigwigs = expand( rules.pileup_multimapping.output.bigwig, sample = samples_list ),
        chromsizes = get_chromsizes,
    output:
        bedgraph = temp( 'analysis/all/multimap-pileup.bedgraph' ),
        bigwig = 'analysis/all/multimap-pileup.bw'
    conda: 'envs/bigwigmerge.yaml'
    resources:
        mem_mb = 8000,
        runtime = '1h',
    params:
        sign=''
    shell:
        '''
        bigWigMerge {input.bigwigs} {output.bedgraph} && \
        awk -v OFS=\"\\t\" '{{print $1,$2,$3, {params.sign}log($4+0.0000000000000000001)/log(10)}}' {output.bedgraph} | \
            LC_COLLATE=C sort -k1,1 -k2,2n -T $TMPDIR -S 8G > {output.bedgraph}.log && \
        bedGraphToBigWig {output.bedgraph}.log {input.chromsizes} {output.bigwig} && \
        rm {output.bedgraph}.log
        '''


use rule aggregate_multimapping as aggregate_primaries with:
    input:
        bigwigs = expand( rules.pileup_primaries.output.bigwig, sample = samples_list ),
        chromsizes = get_chromsizes,
    output:
        bedgraph = temp( 'analysis/all/primary-pileup.bedgraph' ),
        bigwig = 'analysis/all/primary-pileup.bw'
    params:
        sign='-'


rule multimap_ratio_pileup:
    input:
        primary = rules.aggregate_primaries.output.bigwig,
        multimaps = rules.aggregate_multimapping.output.bigwig,
        chromsizes = get_chromsizes,
    output:
        bedgraph = temp('analysis/all/logratio-pileup.bedgraph'),
        bigwig = 'analysis/all/logratio-pileup.bw'
    resources:
        mem_mb = 8000,
        runtime = '1h',
    conda: 'envs/bigwigmerge.yaml'
    shell:
        '''
        bigWigMerge -threshold="-10000" {input.primary} {input.multimaps} {output.bedgraph}.tmp &&\
        LC_COLLATE=C sort -k1,1 -k2,2n -T $TMPDIR -S 8G {output.bedgraph}.tmp > {output.bedgraph} && \
        bedGraphToBigWig {output.bedgraph} {input.chromsizes} {output.bigwig} && \
        rm {output.bedgraph}.tmp
        '''
        

rule multimapping_coverage_stats:
    input:
        bedgraph = rules.multimap_ratio_pileup.output.bedgraph,
        contigs = get_contigs,
        chromsizes = get_chromsizes,
    output:
        'analysis/all/multimap-stats.tsv'
    resources:
        mem_mb = 4096,
        runtime = '1h',
    conda:
        'envs/pandas.yaml'
    params:
        scripts = config['_external_scripts']
    shell:
        'python {params.scripts}/summarize_multimapping.py '
        '--contigs-file {input.contigs} --chromsizes-file {input.chromsizes} '
        '-bg {input.bedgraph} > {output}'
