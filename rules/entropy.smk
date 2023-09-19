

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
            {input.bam} > {output} 2> {log} && \
        samtools faidx {output}
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
        window_size = config['cnv_window_size'],
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
    output:
        'analysis/entropy/regions/{region}.txt'
    params:
        scripts = config['_external_scripts'],
    conda: "envs/bedtools.yaml"
    group: 'alleles'
    params:
        window_size = int( config['cnv_window_size'] )
    shell:
        """
        region_to_bed=$( echo "{wildcards.region}" | awk -F"[:-]" '{{ print $1\"\\t\"$2+1\"\\t\"$3 }}' );
        regions=$(bedtools makewindows -b <( echo $region_to_bed ) -w {params.window_size} | awk '{{ print $1":"$2"-"$3 }}');
        \
        for region in regions; do
            for cons in {input.consensuses}; do samtools faidx -n0 $cons $region | grep -v "^>"; done | \
                python {params.scripts}/allele-stats.py - | \
                awk -v region=$region -v OFS="\t" '{{print region,$0}}' > {output}
        done
        """
        

checkpoint entropy_chunk_genome:
    input:
        get_chromsizes,
    output:
        temp('analysis/entropy/windows.bed')
    conda:
        "envs/bedtools.yaml"
    params:
        window_size = int( config['cnv_window_size']*1000 )
    shell:
        'bedtools makewindows -g {input} -w {params.window_size} > {output}'


def get_chunked_theta_estimates(wildcards):

    regions_file = checkpoints.entropy_chunk_genome.get(**wildcards).output[0]

    regions = []
    with open(regions_file) as f:
        for line in f:
            chrom,start,end = line.strip().split('\t')
            regions.append(
                f'{chrom}:{start}-{end}'
            )

    return expand(
        rules.get_alleles.output, region = regions
    )


rule aggregate_thetas:
    input:
        get_chunked_theta_estimates,
    output:
        'analysis/entropy/theta_estimates.tsv'
    shell:
        'cat {input} > {output}'



