from math import sqrt, log10

## Section 1 - species abundances
#  Count the number of mapped reads falling on each "gene" feature in the gff
#  Save genes x counts for each genome in a separate file
##
rule feature_counts:
    input:
       bam = rules.markduplicates.output,
       bamindex = rules.bam_index.output,
       gff = 'genomes/{genome}/genomic.gff',
    output:
        'analysis/samples/{sample}/feature_counts/{genome}.tsv'
    resources:
        mem_mb = double_on_failure(config['resources']['htseq_count']['mem_mb']),
        runtime = double_on_failure_time(config['resources']['htseq_count']['runtime'])
    threads: 1
    conda:
        'envs/subread.yaml'
    log:
        'logs/feature_count/{sample}_{genome}.log'
    benchmark:
        'benchmark/feature_count/{sample}_{genome}.tsv'
    params:
        scripts = config['_external_scripts'],
        quality = config['min_count_quality']
    group:
        "{sample}_counting"
    shell:
        """
        python {params.scripts}/query-gff -i {input.gff} -type CDS -attr gene -gff > {input.gff}.tmp && \
        featureCounts -a {input.gff}.tmp -o {output} -t CDS -g gene \
            -Q {params.quality} --primary --ignoreDup -p {input.bam} --extraAttributes ID > {log} 2>&1 &&
        rm {input.gff}.tmp
        """

# Aggregate the counts for each genome into a single file 
# which describes relative abundances for each species in the sample
#
rule summarize_abundances:
    input:
        gene_counts = lambda w :expand(rules.feature_counts.output, genome = genomes_list, sample = w.sample)
    output:
        'analysis/samples/{sample}/feature_counts.tsv'
    params:
        genomes = genomes_list
    run:
        import pandas as pd

        def format_featurecounts_output(filename, genome):
            df = pd.read_csv(filename, sep = '\t', skiprows = 1).iloc[:,[0,-1]]\
                .drop_duplicates(subset = ['Geneid'])

            df = df.rename(columns = {'Geneid' : 'gene', df.columns[-1] : genome})\
                .set_index('gene')
            return df

        pd.concat([
                format_featurecounts_output(filename, genome)
                for genome, filename in zip(params.genomes, input.gene_counts)
            ], axis = 1,
        ).fillna(0.).to_csv(output[0], sep = '\t')



rule bamcoverage:
    input:
        bam = rules.filter_bamfile.output,
        bamindex = rules.bam_index.output,
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
        quality = config['min_count_quality']
    shell:
        """
        bedtools genomecov -ibam - -bga < {input.bam} | sort -k1,1 -k2,2n > {output.bedgraph} 2> {log} && \
        bedGraphToBigWig {output.bedgraph} {input.chromsizes} {output.bigwig}
        """


## Section 2 - CNV calling
# Using information from the origin of replication and coverage,
# Call CNVs and PTRs for each sample
##
cnv_outprefix = 'analysis/samples/{sample}/MetaCNV'
rule call_cnvs_and_ptrs:
    input:
        fasta = get_reference,
        bigwig = rules.bamcoverage.output.bigwig,
        #ori = rules.summarize_genomes.output.oris,
    output:
        cnv_info = cnv_outprefix + '.ploidy.bgmtx.gz',
        cnv_calls = cnv_outprefix + '.ploidy.bedgraph.gz',
        ptr = cnv_outprefix + '.PTR.tsv',
        cnv_deviations = cnv_outprefix + '.deviations.bed',
    resources:
        mem_mb = double_on_failure(config['resources']['call_cnvs_and_ptrs']['mem_mb']),
        runtime = double_on_failure_time(config['resources']['call_cnvs_and_ptrs']['runtime'])
    threads: 1
    conda:
        "envs/metaCNV.yaml"
    log:
        'logs/metaCNV/{sample}.log'
    benchmark:
        'benchmark/metaCNV/{sample}.tsv'
    params:
        outprefix = cnv_outprefix,
        scripts = config['_external_scripts']
    shell:
        """
        python {params.scripts}/metaCNV/RunMetaCNV \
            --fasta {input.fasta} \
            --bigwig {input.bigwig} \
            --ori {input.ori} \
            --outprefix {params.outprefix} > {log} 2>&1 && \
        \
        gunzip -c {output.cnv_calls} | \
            awk -v OFS=\"\t\" '$4!=1 {{print $1, $2, $3,\"{wildcards.sample}\",$4}}' > {output.cnv_deviations}
        """



rule aggregate_ploidy_changes:
    input:
        cnv_calls = expand(rules.call_cnvs_and_ptrs.output.cnv_deviations,
                        sample = samples_list),
    output:
        temp('processing/variants/cnv_calls.bed')
    params:
        ploidy = config['base_ploidy']
    shell:
        "cat {input.cnv_calls} | "
        "awk -v OFS=\"\t\" '{{print $1, $2, $3, $4, $5*{params.ploidy} }}' "
        "> {output}"


wildcard_constraints:
    region='[a-zA-Z0-9._-]+:[0-9]+-[0-9]+'

rule callvariants:
    input:
        samples = expand(rules.markduplicates.output, sample = samples_list),
        reference = get_reference,
        cnv_calls = rules.aggregate_ploidy_changes.output,
    output:
        vcf = temp('processing/variants/{region}.vcf'),
    resources:
        mem_mb = double_on_failure(config['resources']['callvariants']['mem_mb']),
        runtime = double_on_failure_time(config['resources']['callvariants']['runtime'])
    threads: config['resources']['callvariants']['threads']
    log:
        'logs/freebayes/{region}.log'
    benchmark:
        'benchmark/freebayes/{region}.tsv'
    params:
        ploidy = config['base_ploidy'],
    conda:
        'envs/freebayes.yaml'
    shell:
        """
        freebayes -f {input.reference} \
            --region {wildcards.region} \
            --ploidy {params.ploidy} \
            --cnv-map {input.cnv_calls} \
            --pooled-discrete \
            --pooled-continuous \
            --use-best-n-alleles 6 \
            --min-alternate-count 2 \
            --min-alternate-fraction 0.01 \
            --haplotype-length 3 \
            --allele-balance-priors-off \
            --report-genotype-likelihood-max \
            --theta 0.01 \
            {input.samples} > {output.vcf}
        """


rule make_snpEff_db:
    input:
        fasta = get_reference,
        gff = get_gff,
    output:
        config = 'analysis/snpEff/snpEff.config',
        done = touch('analysis/snpEff/done.txt')
    conda:
        "envs/snpEff.yaml"
    params:
        db = 'IBD',
        genome_dir = 'analysis/snpEff/data/'
    log:
        'logs/make_snpEff_db.log'
    shell:
        """
        mkdir -p {params.genome_dir}/{params.db} && \
        echo -e "{params.db}.genome : {params.db}" > {output.config} && \
        cp {input.gff} {params.genome_dir}/{params.db}/genes.gff && \
        cp {input.fasta} {params.genome_dir}/{params.db}/sequences.fa && \
        snpEff build -gff3 -v -c {output.config} -noCheckCDS -noCheckProtein {params.db} > {log} 2>&1
        """


rule annotate_vcf:
    input:
        vcf = rules.callvariants.output,
        snpEff_database = rules.make_snpEff_db.output,
        snpEff_config = rules.make_snpEff_db.output.config,
    output:
        temp('processing/annotation/{region}.vcf')
    conda:
        "envs/snpEff.yaml"
    log:
        'logs/annotate/{region}.log'
    benchmark:
        'benchmark/annotate/{region}.tsv'
    resources:
        mem_mb = double_on_failure(config['resources']['annotate_vcf']['mem_mb']),
        runtime = double_on_failure_time(config['resources']['annotate_vcf']['runtime'])
    threads: config['resources']['annotate_vcf']['threads']
    params:
        config = 'analysis/snpEff/snpEff.config',
    shell:
        "awk '$4 !~ /N/' {input.vcf} | snpEff -c {input.snpEff_config} IBD > {output}"


checkpoint chunk_genome:
    input:
        get_chromsizes,
    output:
        temp('genomes/all/windows.bed')
    conda:
        "envs/bedtools.yaml"
    params:
        window_size = int( config['variant_calling_window_size'] * 100000 )
    shell:
        'bedtools makewindows -g {input} -w {params.window_size} > {output}'


def get_chunked_variant_calls(wildcards):

    regions_file = checkpoints.chunk_genome.get(**wildcards).output[0]

    with open(regions_file) as f:
        regions = [
            '{}:{}-{}'.format(*map(str, line.strip().split('\t')[:3])) 
            for line in f
        ]

    return expand(
        rules.annotate_vcf.output, region = regions
    ) 


def fdr_to_qual(fdr):
    return int(-10 * log10(fdr))

rule merge_vcfs:
    input:
        vcfs = get_chunked_variant_calls,
        reference = get_reference
    output:
        vcf = protected('analysis/all/variants.bcf')
    conda:
        'envs/bcftools.yaml'
    params:
        min_reads = config['min_supporting_reads'],
        qual = fdr_to_qual(config['variant_fdr']),
    conda:
        "envs/bcftools.yaml"
    log:
        "logs/merge_vcfs.log"
    resources:
        mem_mb = double_on_failure(config['resources']['merge_vcfs']['mem_mb']),
        runtime = double_on_failure_time(config['resources']['merge_vcfs']['runtime'])
    threads: config['resources']['merge_vcfs']['threads']
    shell:
        """
        bcftools concat {input.vcfs} | \
        bcftools norm -d none -f {input.reference} | \
        bcftools filter -i 'QUAL>={params.qual}' -Ov > {output.vcf} 2> {log}
        """


rule get_sample_vcf:
    input:
        rules.merge_vcfs.output
    output:
        temp('analysis/samples/{sample}/vcf.bcf')
    conda:
        "envs/bcftools.yaml"
    shell:
        "bcftools view {input} -c1 -Oz -s {wildcards.sample} > {output} && bcftools index {output}"


