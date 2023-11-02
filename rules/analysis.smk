from math import sqrt, log10

## Section 1 - species abundances
#  Count the number of mapped reads falling on each "gene" feature in the gff
#  Save genes x counts for each genome in a separate file
##
rule get_CDS_features:
    input:
        'genomes/{genome}/genomic.gff'
    output:
        temp('genomes/{genome}/cds.gff')
    log:
        'logs/{genome}/CDS_features.log'
    params:
        scripts = config['_external_scripts']
    shell:
        'python {params.scripts}/query-gff -i {input} -type CDS -attr gene -gff > {output}'


rule feature_counts:
    input:
       bam = get_bam,
       bamindex = get_bai,
       gff = rules.get_CDS_features.output,
    output:
        temp('analysis/samples/{sample}/feature_counts/{genome}.tsv')
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
        quality = config['min_count_quality'],
        paired_param = lambda wildcards : '-p' if config['samples'][wildcards.sample]['is_paired'] else '',
    group:
        "{sample}_counting"
    shell:
        """
        featureCounts -a {input.gff} -o {output} -t CDS -g gene \
            -Q {params.quality} --primary --ignoreDup {params.paired_param} {input.bam} --extraAttributes ID > {log} 2>&1
        """

# Aggregate the counts for each genome into a single file 
# which describes relative abundances for each species in the sample
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


##
# Section 2: Variant calling
# These rules split the genomes into large chunks, and call variants across all samples 
# on these large chunks. The VCFs from each chunk are filtered for low-confidence variants
# and annotated using SNPEff. Finally, the VCFs are merged into the "master" variant call list.
##

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


wildcard_constraints:
    region='[a-zA-Z0-9._-]+:[0-9]+-[0-9]+'


rule callvariants:
    input:
        samples = expand(rules.markduplicates.output.bam, sample = samples_list),
        reference = get_reference,
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
            --region "{wildcards.region}" \
            --ploidy {params.ploidy} \
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

##
# Split the genome into chunks for variant calling
##
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

##
# Run the rule above, and read the chunks.
# Then, return a list of target VCFs, one per region to
# instruct the pipeline to generate each.
##
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

##
# This rule takes each VCF chunk, merges them, and filters out low-quality variants.
## 
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


