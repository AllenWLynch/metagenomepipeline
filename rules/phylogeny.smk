
import os


rule get_alignment_region:
    input:
        rules.move_genome.output.gff,
    output:
        regions = temp('analysis/MSA/{genome}/{core_gene}.bed'),
    params:
        scripts = config['_external_scripts'],
    group: 'consensus'
    shell:
        # The first step sorts the contigs file by chromosome, then selects contigs for the given genome
        # Next, the gff file is joined to the contigs file,
        # and then the gff file is filtered to only include core genes which is outputted as a bedfile with an informative name
        """
        python {params.scripts}/query-gff \
            -i {input}
            -type gene \
            -attr gene \
            -val {wildcards.core_gene} \
            -f '{chr}\t{start}\t{end}\tchr={chr};start={start};end={end};genome={wildcards.genome};gene={attributes[gene]};id={attributes[ID]}\n' | \
            head -n1 > {output}
        """


def format_sample_metadata(sampledata):

    def flat_dict(dict, sep = ';', rel = '='):
        return sep.join([rel.join([k, str(v)]) for k,v in dict.items()])

    metadata = sampledata.pop('metadata')
    sampledata.update(
        flat_dict({'metadata.' + k : v for k,v in metadata.items()})
    )

    return flat_dict(sampledata)


rule get_dominant_strain_genotype:
    input:
        vcf = get_vcf_input,
        regions = rules.get_alignment_region.output.regions,
        reference = rules.move_genome.output.fasta,
    output:
        vcf = temp('analysis/samples/{sample}/consensus/{genome}/{core_gene}_variants.bcf'),
        consensus = temp('analysis/sample/{sample}/consensus/{genome}/{core_gene}_consensus.fna'),
    conda:
        'envs/consensus.yaml'
    group: 'consensus'
    params:
        metadata = lambda w : format_sample_metadata(config['samples'][w.sample]),
    shell:
        """
        bcftools filter --IndelGap 5 {input.vcf} --regions-file {input.regions} | \
        bcftools norm -m- | \
        bcftools view -i 'FORMAT/AO>(0.5*INFO/DP)' -Oz > {output.vcf} \
        && bcftools index {output.vcf} && \
        \
        awk -v OFS=\"\" '{{print $0,\";sample={wildcards.sample};{params.metadata}\"}}' | \
        bedtools get-fasta -f {input.reference} -bed - -name | \
        bcftools consensus -s - -f - {output.vcf} > {output.consensus}
        """


rule get_reference_sequence:
    input:
        reference = rules.move_genome.output.fasta
        regions = rules.get_alignment_regions.output.regions
    output:
        'genomes/{genome}/core_genes/{core_gene}.fna',
    group: 'consensus'
    shell:
        'bedtools get-fasta -f {input.reference} -bed {input.regions} -name > {output.consensus}'


## 
# In this step, load in the abundance matrices for core genes for each sample, and then
# choose a list of samples to include in the phylogeny. This is done by choosing
# a representative sample per subject, and by filtering out samples with too low coverage.
## 
checkpoint get_most_abundant_samples:
    input:
        lambda w : expand(rules.summarize_abundances.output, sample =samples_list),
    output:
        'analysis/MSA/{genome}/samples_list.txt'
    shell:
        pass


def get_usable_samples(wildcards):

    checkpoint_output = checkpoints.get_most_abundant_samples.get(**wildcards).output[0]
    
    with open(checkpoint_output, 'r') as f:
        sample_ids = [x.strip() for x in f]

    return expand(
        rules.get_sample_vcf.output, sample = sample_ids
    )

##
# A tad messy, but here, "genome" refers to the reference genome for a particular species,
# while "alternate_ref" refers to other example genomes for that species.
# In the rule below, we must trade off "genome" for "alternate_ref" in the input names.
## 
rule make_msa:
    input:
        sample_consensus = get_usable_samples,
        reference_sequences = lambda wildcards : expand(rules.get_reference_sequence.output, 
                                    genome = config['sequences'][wildcards.genome]['alternate_refs']) + get_genome_reference(wildcards.genome),
    output:
        msa = 'analysis/MSA/{genome}/{core_gene}.msa.fna',
    shell:
        "echo HERE"



rule make_tree:
    input:
        expand(rules.make_msa.output.msa, core_gene=config['core_genes']),
    output:
        tree = 'analysis/MSA/{genome}/tree.nwk',
    shell:
        "echo HERE"
