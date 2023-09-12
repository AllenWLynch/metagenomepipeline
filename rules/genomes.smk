
wildcard_constraints:
    genome = '[a-zA-Z0-9._-]+|all'

# Fetch the genome from the same NCBI endpoint.
download_url = '"https://api.ncbi.nlm.nih.gov/datasets/v2alpha/genome/accession/{genome}/download?include_annotation_type={filetype},SEQUENCE_REPORT&filename={genome}_fasta.zip"'
rule download_genome:
    output:
        temp('{genome}_fasta.zip')
    params:
        url = lambda wildcards : download_url.format(genome = wildcards.genome, filetype='GENOME_FASTA')
    group: 'download'
    shell:
        'curl -OJX GET {params.url} --fail -H "Accept: application/zip"'


rule download_fasta:
    output:
        'genomes/{genome}/genomic.fa'
    group: 'download'
    params:
        url = lambda wildcards : 'ftp://ftp.ebi.ac.uk/pub/databases/ena/wgs/public/cab/{ENA_ID}.fasta.gz'\
                                    .format(ENA_ID = config['genomes'][wildcards.genome]['ENA_ID'])
    shell:
        'curl -OJX GET {params.url} --fail && -o {output}.gz && gunzip {output}.gz'


rule download_gff:
    output:
        'genomes/{genome}/genomic.gff'
    group: 'download'
    params:
        scripts = config['_external_scripts'],
        genome_id = lambda w : config['genomes'][w.genome]['MGYG'],
        gut_id = lambda w : config['genomes'][w.genome]['GUT_ID'],
    shell:
        'bash {params.scripts}/get-uhgg-gff {params.genome_id} {params.gut_id}.gff.gz && '
        'gunzip -c {params.gut_id}.gff.gz > {output} && '
        'rm {params.gut_id}.gff.gz'


# The NCBI endpoint returns a zip file with a nasty directory scheme.
# This rule normalizes the directory structure and moves the genome and GFF to the top level.
# Then, it extracts and normalizes the genome fasta.
rule unpack_genome:
    input:
        fasta_download = rules.download_genome.output,
        gff_download = rules.download_gff.output,
    output:
        fasta = 'genomes/{genome}/genomic.fa',
        gff = 'genomes/{genome}/genomic.gff',
        assembly_data_report = 'genomes/{genome}/assembly_data_report.jsonl'
    group: 'download'
    params:
        fasta_download = lambda wildcards : \
            'genomes/{genome}/{genome}_{name}_genomic.fna'.format(genome = wildcards.genome, name = "*"),
    shell:
        '''
        mkdir -p genomes/tmp/{wildcards.genome} && mkdir -p genomes/{wildcards.genome} && \
        unzip -o {input.fasta_download} -d genomes/tmp/{wildcards.genome} && \
        mv genomes/tmp/{wildcards.genome}/ncbi_dataset/data/{wildcards.genome}/* genomes/{wildcards.genome}/ && \
        mv {params.fasta_download} {output.fasta} && \
        mv genomes/tmp/{wildcards.genome}/ncbi_dataset/data/assembly_data_report.jsonl genomes/{wildcards.genome}/ && \
        rm -rf genomes/tmp/{wildcards.genome} && \
        \
        gunzip -c {input.gff_download} > genomes/{wildcards.genome}/genomic.gff
        '''


# This rule summarizes the merged genome for use in downstream analysis.
# Some tools require a fasta index, chromosome sizes, or a genome dictionary.
rule summarize_genome:
    input:
        'genomes/{genome}/genomic.fa',
    output:
        index = 'genomes/{genome}/genomic.fa.fai',
        chromsizes = 'genomes/{genome}/genomic.chromsizes',
        dict = 'genomes/{genome}/genomic.dict',
    conda:
        "envs/samtools.yaml"
    group: "download"
    shell:
        '''
        samtools faidx {input} && \
        cut -f1,2 {input}.fai > {output.chromsizes} && \
        samtools dict {input} -o {output.dict}
        '''


rule make_contigs_file:
    input:
        fasta = rules.unpack_genome.output.fasta,
        assembly_data_report = rules.unpack_genome.output.assembly_data_report,
    output:
        'genomes/{genome}/contigs.tsv'
    group:
        'download'
    params:
        header = '#' + '\t#'.join('genome, contig, species'.split(', ')),
        species = lambda wildcards : config['genomes'][wildcards.genome]['species'],
    shell:
        'echo -e "{params.header}" > {output} && ' 
        'grep ">" {input.fasta} | cut -c2- | cut -f1 -d" " | awk -v OFS="\t" \'{{print "{wildcards.genome}", $0, "{params.species}"}}\' >> {output}'
        

'''
grep ">" {input} | cut -c2- | tr " " "\\t" | awk -v OFS="\\t" \'{{print "{wildcards.genome}", $0}}\' | cut -f1-6 > {output.contigs}.tmp1 && \
grep ">" {input} | cut -f6- -d" " > {output.contigs}.tmp2 && \
\
echo -e \"{params.header}\" > {output.contigs} && \
paste {output.contigs}.tmp1 {output.contigs}.tmp2 >> {output.contigs} && \
rm {output.contigs}.tmp1 {output.contigs}.tmp2
'''
#gc_skew_bedgraph = 'genomes/gc_skew.bedgraph.gz',
 #       oris = 'genomes/replication_factors.bed'
"""bash {params.scripts}/ORIfinder \
            {input} {output.chromsizes} {output.gc_skew_bedgraph} > {output.oris}"""