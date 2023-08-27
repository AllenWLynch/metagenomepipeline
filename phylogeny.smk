

rule get_alignment_regions:
    output:
        regions = 'analysis/MSA/regions.txt',
        orthogroups = 'analysis/MSA/orthogrous.txt'
    run:
        pass


rule get_dominant_strain_genotype:
    input:
        vcf = rules.get_sample_vcf.output,
        reference = get_reference,
        regions = rules.get_alignment_regions.output.regions,
    output:
        vcf = temp('analysis/samples/{sample}/dominant_strain_variants.bcf'),
        consensus = 'analysis/samples/{sample}/consensus.fna',
    shell:
        """
        bcftools filter --IndelGap 5 {input} --regions {input.regions} | \
        bcftools norm -m- | \
        bcftools view -i 'AF>0.5' -Oz > {output.vcf} \
        && bcftools index {output.vcf} && \
        samtools faidx <(cat {input.regions}) {params.regions_list} | \
        bcftools consensus -s - -f - {output.vcf} > {output.consensus}
        """