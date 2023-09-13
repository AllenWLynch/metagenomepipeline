
wildcard_constraints:
    genome = '[a-zA-Z0-9._-]+|all'


rule download_annotation:
    output:
        temp('genomes/{genome}/annotation.gff')
    group: 'download'
    params:
        scripts = config['_external_scripts'],
        genome_id = lambda w : config['genomes'][w.genome]['MGYG'],
        gut_id = lambda w : config['genomes'][w.genome]['GUT_ID'],
    shell:
        'bash {params.scripts}/get-uhgg-gff {params.genome_id} {params.gut_id}.gff.gz && '
        'gunzip {params.gut_id}.gff.gz -c > {output} && '
        'rm {params.gut_id}.gff.gz'


rule genome:
    input:
        rules.download_annotation.output
    output:
        fasta = 'genomes/{genome}/genomic.fa',
        gff = 'genomes/{genome}/genomic.gff',
    run:
        in_fasta = False
        with open(input[0], 'r') as f, open(output.fasta, 'w') as fasta, open(output.gff, 'w') as gff:
            for line in f:
                if line.startswith('##FASTA'):
                    in_fasta = True
                elif in_fasta:
                    fasta.write(line)
                else:
                    gff.write(line)


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
        'genomes/{genome}/genomic.fa'
    output:
        'genomes/{genome}/contigs.tsv'
    group:
        'download'
    params:
        header = '#' + '\t#'.join('genome, contig, species'.split(', ')),
        species = lambda wildcards : config['genomes'][wildcards.genome]['species'].replace(' ','-'),
    shell:
        'echo -e "{params.header}" > {output} && ' 
        'grep ">" {input} | cut -c2- | cut -f1 -d" " | awk -v OFS="\t" \'{{print "{wildcards.genome}", $0, "{params.species}"}}\' >> {output}'
        

#gc_skew_bedgraph = 'genomes/gc_skew.bedgraph.gz',
 #       oris = 'genomes/replication_factors.bed'
"""bash {params.scripts}/ORIfinder \
            {input} {output.chromsizes} {output.gc_skew_bedgraph} > {output.oris}"""