

rule get_phylogeny_gff:
    input:  
        get_gff
    output:
        temp('analysis/phylogeny/tmp/gff.gff')
    params:
        genes = config['phylogeny_genes']
    shell:
        "echo HERE"


rule get_consensus:
    input:
        rules.get_sample_vcf.output
    output:
        protected('analysis/samples/{sample}/consensus/{element}.fna')
    params:
        region = get_region,
    conda:
        "envs/bcftools.yaml"
    shell:
        """
        bcftools norm -m+any {input} --region {params.region} | \
        bcftools query -f '%CHROM %POS %REF %ALT %INFO/AC %INFO/DP' | \
        python scripts/consensus.py --region {params.region} > {output}
        """

rule get_consensus_sequence:
    input:
        rules.get_allele_counts.output
    output:
        protected('analysis/samples/{sample}/consensus/{element}.fna')
    params:
        region = get_region,
    shell:
        """
        python scripts/consensus.py {input} \
            -ref {input.reference} \
            --region {params.region} > {output}
        """


rule get_phylogeny_fasta:
    input:
        gff = rules.get_phylogeny_gff.output
        consensus_genomes = expand(
            rules.get_consensus.output,
            sample = config['samples'].keys()
        )
    output:
        temp('analysis/phylogeny/tmp/fasta.fa')
    shell:
        """
        for consensus in ({input.consensus_genomes}); do \
            bedtools getfasta -f $consensus -g {input.gff} --name > analysis/phylogeny/tmp/{$(basename $consensus)}; \
        endfor && \
        cat {input.consensus_genomes} > {output}
        """

rule get_phylogeny_OGs:
    output: temp('analysis/phylogeny/tmp/OGs.tsv')


rule make_phylogeny:
    input:
        OGs = rules.get_phylogeny_OGs.output,
        fasta = rules.get_phylogeny_fasta.output,
    output:
        protected('analysis/phylogeny/tmp/fasta.fa')


"""
python scripts/get_variant_array.py \
    -vcf analysis/samples/{sample}/vcf.bcf \
    -ref analysis/samples/{sample}/consensus.fna \

"""

'''rule summarize_small_variants_genome:
    input:
        vcf = rules.get_sample_vcf.output,
        coverage = rules.bamcoverage.output,
        reference = get_genome_reference,
        gff = get_genome_gff,
    output:
        protected(directory('genomes/{genome}/{sample}/allele_counts/'))
    params:
        output_dir = lambda w : summarize_variant_path.format(**w)
    shell:
        """
        python scripts/get_consensus.py \
            -vcf {input.vcf} \
            -coverage {input.coverage} \
            -ref {input.reference} \
            -gff {input.gff} \
            -o {params.output_dir}
        """

rule summarize_small_variants_gene:
    input:
        rules.summarize_small_variants_genome.output
    output:
        protected('gene_arrays/{sample}-{species}-{gene}-variant_counts.mm')
    params:
        gene_start, gene_end, strand
    conda:
        "envs/scipy.yaml"
    shell:
        """
        python -c 'import scipy; mtx = scipy.io.mmread('{input}').tocsc(); \
                mtx = mtx[:, {params.gene_start}:{params.gene_end}:(2*{params.strand} - 1) ]; \
                mtx.eliminate_zeros(); scipy.io.mmwrite('{output}', mtx)'
        """


def list_species():
    return config['genomes'].keys()

rule generate_phylogenetic_tree:
    input:
        expand(rules.get_consensus_strain_genes, 
                sample = config['samples'].keys(), 
                species = list_species, 
                gene = config['phylogeny_genes']
        )
    output:
        protected('MSA/MSA.tree')
    shell:
        "echo here"

'''