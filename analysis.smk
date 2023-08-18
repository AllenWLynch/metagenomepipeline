
rule htseq_count:
    input:
       bam = rules.markduplicates.output,
       gff = get_gff
    output:
        protected('analysis/samples/{sample}/htseq.counts')
    resources:
        mem_mb = double_on_failure(config['resources']['htseq_count']['mem_mb']),
        runtime = double_on_failure(config['resources']['htseq_count']['runtime'])
    threads: 1
    conda:
        'envs/htseq.yaml'
    shell:
        'htseq-count -i Name -r pos -s no -a 0 --secondary-alignments ignore -t gene {input.bam} {input.gff} > {output}'


rule bamcoverage:
    input:
        bam = rules.markduplicates.output,
        chromsizes = rules.summarize_genomes.output.chromsizes
    output:
        bedgraph = temp('analysis/samples/{sample}/coverage.bedgraph'),
        bigwig = protected('analysis/samples/{sample}/coverage.bigwig')
    log:
        'logs/coverage/{sample}.coverage.log',
    resources:
        mem_mb = double_on_failure(config['resources']['bamcoverage']['mem_mb']),
        runtime = double_on_failure(config['resources']['bamcoverage']['runtime'])
    threads: config['resources']['bamcoverage']['threads']
    conda:
        'envs/bedtools.yaml'
    shell:
        """
        samtools view -q 0 -b -F 0x400 -F 0x100 -F 0x800 {input.bam} | \
        bedtools genomecov -ibam - -bga > {output.bedgraph} && \
        bedGraphToBigWig {output.bedgraph} {input.chromsizes} {output.bigwig}
        """


def list_input_samples(w):
    return [rules.markduplicates.output[0].format(sample = sample)
            for sample in config['groups'][w.group]
        ]

rule callvariants:
    input:
        samples = list_input_samples,
        reference = get_reference,
    output:
        vcf = 'analysis/groups/{group}/vcf.gz',
        #stats = protected('analysis/groups/{group}/vcf.stats')
    log:
        'logs/mutect/{group}.mutect.log'
    resources:
        mem_mb = double_on_failure(config['resources']['callvariants']['mem_mb']),
        runtime = double_on_failure(config['resources']['callvariants']['runtime'])
    threads: config['resources']['callvariants']['threads']
    params:
        input_list = lambda w : '-I ' + ' -I '.join(list_input_samples(w)),
        imputed_af = 1/3, #len(config['groups']),
        gatk_path = config['gatk_path']
    #conda:
    #    'envs/gatk.yaml'
    shell:
        """
        {params.gatk_path} \
        Mutect2 -O {output.vcf} \
            {params.input_list} \
            --reference {input.reference} \
            --af-of-alleles-not-in-resource {params.imputed_af}
        """

sample_group_dict = {sample : group for group, samples in config['groups'].items()
                     for sample in samples}

# splits samples and normalizes VCF to remove multiallelic sites
rule get_sample_vcf:
    input:
        vcf = lambda w : rules.callvariants.output.vcf.\
            format(group = sample_group_dict[w.sample]),
        reference = get_reference
    output:
        protected('analysis/samples/{sample}/vcf.bcf')
    conda:
        "envs/bcftools.yaml"
    shell:
        "gunzip -c {input.vcf} | bcftools view -c1 -Ov -s {wildcards.sample} | "
        "bcftools norm -f {input.reference} -Oz > {output}"


rule get_allele_counts:
    input:
        rules.get_sample_vcf.output
    output:
        protected('analysis/samples/{sample}/allele_counts.tsv.gz')
    conda:
        "envs/bcftools.yaml"
    shell:
        """
        bcftools norm -m+any -a --exclude 'TYPE=\"indel\"' {input} | \
        bcftools query -f '%CHROM %POS %REF %ALT %INFO/AC %INFO/DP' | \
        python scripts/get_readsupport.py | bgzip -c > {output} && \
        tabix -s1 -b2 -e2 {output}
        """


rule make_snpEff_db:
    input:
        reference = get_reference,
        gff = get_gff,
    output:
        'analysis/snpEff/snpEff.config'
    #conda:
    #    "envs/snpEff.yaml"
    shell:
        "snpEff build -gff3 -v {input.reference} > {output}"


rule annotate_vcf:
    input:
        vcf = rules.get_sample_vcf.output,
        snpEff_database = rules.make_snpEff_db.output,
    output:
        protected('analysis/samples/{sample}/vcf.annotated.bcf')
    #conda:
    #    "envs/snpEff.yaml"
    shell:
        "snpEff "


rule merge_vcfs:
    input:
        expand(rules.annotate_vcf.output, sample = config['samples'].keys())
    output:
        protected('analysis/all_variants.bcf')
    conda:
        "envs/bcftools.yaml"
    shell:
        "bcftools merge -Oz {input} -m+ > {output}"
