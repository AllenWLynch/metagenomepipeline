

rule get_consensus:
    input:
        bam = get_bam,
    output:
        'analysis/samples/{sample}/consensus.fa'
    log:
        'logs/consensus/{sample}.log'
    conda:
        "envs/samtools.yaml"
    shell:
        '''
        samtools consensus \
            -f fasta \
            --mode simple \
            -aa --show-ins no --show-del yes \
            --use-qual \
            --min-MQ 20 \
            --min-depth 5 \
            --call-frac 0.6 \
            {input.bam} > {output} 2> {log}
        '''


rule get_entropy:
    input:
        bam = get_bam,
        fasta = get_reference,
    output:
        'analysis/samples/{sample}/entropy.bed'
    log:
        'logs/entropy/{sample}.log'
    benchmark:
        'benchmark/entropy/{sample}.tsv'
    params:
        scripts = config['_external_scripts'],
        window_size = 200,
        min_depth = 5,
    conda:
        "envs/samtools.yaml"
    shell:
        '''
        bash {params.scripts}/empirical_entropy \
            {input.bam} {input.fasta} {params.window_size} {params.min_depth} | \
            bgzip > {output} 2> {output}.log && \
        tabix -p bed {output}
        '''


rule get_alleles:
    input:
        consensuses = expand(rules.get_consensus.output, sample = samples_list),
        chromsizes = get_chromsizes,
    output:
        'analysis/entropy/alleles.tsv'
    params:
        scripts = config['_external_scripts'],
        window_size = 200,
    shell:
        '''
        echo -e "region\ttheta_mle\tn_alleles\tn_samples\tconsensus_allele\tn_variants_distribution" > {output};
        regions=$(bedtools makewindows -g {input.chromsizes} -w {params.window_size} | awk '{{print $1":"$2 + 1"-"$3}}');
        for region in $regions; do
            for consensus in {input.consensuses}; do samtools faidx $consensus $region | grep -v "^>"; done | \
            python {params.scripts}/allele-stats.py | \
            awk -v region=$region -v OFS="\t" '{{print region,$0}}' \
            >> {output};
        done
        '''

