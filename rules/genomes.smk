

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
rule move_genome:
    input:
        rules.download_genome.output
    output:
        _dir = directory('genomes/{genome}'),
        gff = 'genomes/{genome}/genomic.gff',
        fasta = 'genomes/{genome}/genomic.fa',
        contigs = temp('genomes/{genome}/contigs.tsv'),
    group: 'download'
    params:
        fasta_download = lambda wildcards : \
            'genomes/{genome}/{genome}_{name}_genomic.fna'.format(genome = wildcards.genome, name = "*"),
    shell:
        'mkdir -p genomes/tmp/{wildcards.genome} && mkdir -p genomes/{wildcards.genome} && '
        'unzip -o {wildcards.genome}.zip -d genomes/tmp/{wildcards.genome} && '
        'mv genomes/tmp/{wildcards.genome}/ncbi_dataset/data/{wildcards.genome}/* genomes/{wildcards.genome}/ && '
        'mv {params.fasta_download} {output.fasta} && '
        'cat {output.fasta} | grep ">" | cut -c2- | tr " " "\\t" | awk -v OFS="\\t" \'{{print "{wildcards.genome}", $0}}\' | cut -f1-7 > {output.contigs}'


# This rule merges all the GFF annotations into a single file
# for alignment indexing.
rule merge_annotations:
    input:
        fastas = expand(rules.move_genome.output.fasta, genome = genomes_list),
        gff = expand(rules.move_genome.output.gff, genome = genomes_list),
        contigs = expand(rules.move_genome.output.contigs, genome = genomes_list),
    output:
        fasta = 'genomes/fasta.fna',
        gff = 'genomes/gff.gff',
        contigs = 'genomes/contigs.tsv'
    group : "merge-genomes"
    params:
        header = '#' + '\t#'.join('genome, contig, genus, species, source, strain, type'.split(', '))
    shell:
        "cat {input.fastas} > {output.fasta} && " \
        "cat {input.gff} | grep -v '#' | sort -k1,1 -k5,5n > {output.gff} && " \
        "echo -e \"{params.header}\" > {output.contigs} && cat {input.contigs} >> {output.contigs}"


# This rule summarizes the merged genome for use in downstream analysis.
# Some tools require a fasta index, chromosome sizes, or a genome dictionary.
rule summarize_genomes:
    input:
        rules.merge_annotations.output.fasta,
    output:
        index = 'genomes/fasta.fna.fai',
        chromsizes = 'genomes/fasta.chromsizes',
        dict = 'genomes/fasta.dict',
        gc_skew_bedgraph = 'genomes/gc_skew.bedgraph.gz',
        oris = 'genomes/replication_factors.bed'
    conda:
        "envs/bedtools.yaml"
    group: "merge-genomes"
    params:
        scripts = config['_external_scripts']
    shell:
        """
        samtools faidx {input} && \
        cut -f1,2 {output.index} > {output.chromsizes} && \
        samtools dict {input} -o {output.dict} && \
        bash {params.scripts}/ORIfinder \
            {input} {output.chromsizes} {output.gc_skew_bedgraph} > {output.oris}
        """