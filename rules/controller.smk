
import re

# Define conveniance methods for accessing the fasta and gff files for downstream analysis.
# This allows you to provide your own reference or gff instead of automatically downloading from NCBI. 
def get_reference(wildcards):
    return rules.merge_annotations.output.fasta

def get_gff(wildcards):
    return rules.merge_annotations.output.gff

def get_contigs(wildcards):
    return 'genomes/all/contigs.tsv'

def get_chromsizes(wildcards):
    return 'genomes/all/genomic.chromsizes'

def get_bam(wildcards):
    return rules.markduplicates.output.bam


def double_on_failure(base_resources):
    def _double_on_failure(wildcards, attempt):
        return base_resources * 2**(attempt - 1)
    return _double_on_failure


def double_on_failure_time(base_time):

    match = re.match(r'(\d+)([smhd])', base_time)
    if not match:
        raise ValueError("Invalid input format")

    value, unit = match.groups()
    value = int(value)

    def _double_on_failure_time(wildcards, attempt):
        return f'{int(value * 2**(attempt - 1))}{unit}'
    
    return _double_on_failure_time

try:
    genomes_list = list(config['genomes'].keys())
    samples_list = list(config['samples'].keys())
except KeyError:
    genomes_list = []
    samples_list = []


include: "genomes.smk"

if config['_run_pipeline'] in ['build-ref','variants','align']:    
    include: "alignment.smk"
    include: "postprocessing.smk"

if config['_run_pipeline'] == 'build-ref':

    targets = [
        *expand(rules.summarize_genome.output.chromsizes, genome = genomes_list + ['all']),\
        *rules.merge_annotations.output,
        rules.make_index.output,
        rules.make_snpEff_db.output,
    ]

elif config['_run_pipeline'] == 'download-genome':

    targets = [
        expand(rules.summarize_genome.output.chromsizes, genome = config['_download_genomes']),
    ]

elif config['_run_pipeline'] == 'align':

    targets = [
        *expand(rules.markduplicates.output, sample = samples_list),
        *expand(rules.get_multimap_stats.output.stats, sample = samples_list)
    ]

elif config['_run_pipeline'] == 'variants':

    include: "analysis.smk"
    
    targets = [
        rules.merge_vcfs.output,
        *expand(rules.summarize_abundances.output, sample = samples_list),
    ]


rule all:
    input: targets
        
