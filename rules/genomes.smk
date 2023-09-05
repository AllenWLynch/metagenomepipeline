
wildcard_constraints:
    genome = '|'.join(genomes_list + ['all'])


# Fetch the genome from the same NCBI endpoint.
# This keeps genomes and GFF annotations consistent so the pipeline doesn't break unexpectedly.
download_url = '"https://api.ncbi.nlm.nih.gov/datasets/v2alpha/genome/accession/{genome}/download?include_annotation_type=GENOME_FASTA,GENOME_GFF,RNA_FASTA,CDS_FASTA,PROT_FASTA,SEQUENCE_REPORT&filename={genome}.zip"'
rule download_genome:
    output:
        temp('{genome}.zip')
    params:
        url = lambda wildcards : download_url.format(genome = wildcards.genome)
    group: 'download'
    shell:
        'curl -OJX GET {params.url} -H "Accept: application/zip"'


# The NCBI endpoint returns a zip file with a nasty directory scheme.
# This rule normalizes the directory structure and moves the genome and GFF to the top level.
# Then, it extracts and normalizes the genome fasta and creates a contigs file for use in downstream analysis.
rule move_genome:
    input:
        rules.download_genome.output
    output:
        gff = 'genomes/{genome}/genomic.gff',
        fasta = 'genomes/{genome}/genomic.fa',
    group: 'download'
    params:
        fasta_download = lambda wildcards : \
            'genomes/{genome}/{genome}_{name}_genomic.fna'.format(genome = wildcards.genome, name = "*"),
    shell:
        'mkdir -p genomes/tmp/{wildcards.genome} && mkdir -p genomes/{wildcards.genome} && '
        'unzip -o {wildcards.genome}.zip -d genomes/tmp/{wildcards.genome} && '
        'mv genomes/tmp/{wildcards.genome}/ncbi_dataset/data/{wildcards.genome}/* genomes/{wildcards.genome}/ && '
        'mv {params.fasta_download} {output.fasta} '


# This rule summarizes the merged genome for use in downstream analysis.
# Some tools require a fasta index, chromosome sizes, or a genome dictionary.
rule summarize_genome:
    input:
        'genomes/{genome}/genomic.fa',
    output:
        index = 'genomes/{genome}/genomic.fa.fai',
        chromsizes = 'genomes/{genome}/genomic.chromsizes',
        dict = 'genomes/{genome}/genomic.dict',
        contigs = 'genomes/{genome}/contigs.tsv'
    conda:
        "envs/bedtools.yaml"
    group: "download"
    params:
        header = '#' + '\t#'.join('genome, contig, genus, species, source, strain, type'.split(', '))
    shell:
        '''
        samtools faidx {input} && \
        cut -f1,2 {input}.fai > {output.chromsizes} && \
        samtools dict {input} -o {output.dict} && \
        \
        grep ">" {input} | cut -c2- | tr " " "\\t" | awk -v OFS="\\t" \'{{print "{wildcards.genome}", $0}}\' | cut -f1-6 > {output.contigs}.tmp1 && \
        grep ">" {input} | cut -f6- -d" " > {output.contigs}.tmp2 && \
        \
        echo -e \"{params.header}\" > {output.contigs} && \
        paste {output.contigs}.tmp1 {output.contigs}.tmp2 >> {output.contigs} && \
        rm {output.contigs}.tmp1 {output.contigs}.tmp2
        '''


# This rule merges all the GFF annotations into a single file
# for alignment indexing.
rule merge_annotations:
    input:
        fastas = expand(rules.move_genome.output.fasta, genome = genomes_list),
        gff = expand(rules.move_genome.output.gff, genome = genomes_list),
    output:
        fasta = 'genomes/all/genomic.fa',
        gff = 'genomes/all/genomic.gff',
    shell:
        "cat {input.fastas} > {output.fasta} && " \
        "cat {input.gff} | grep -v '#' | sort -k1,1 -k5,5n > {output.gff} && " \


#gc_skew_bedgraph = 'genomes/gc_skew.bedgraph.gz',
 #       oris = 'genomes/replication_factors.bed'
"""bash {params.scripts}/ORIfinder \
            {input} {output.chromsizes} {output.gc_skew_bedgraph} > {output.oris}"""