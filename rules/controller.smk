
import re

genomes_list = [metadata['reference'] for species, metadata in config['genomes'].items()]
samples_list = list(config['samples'].keys())

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

### DELETE THIS !!! ##
# config['groups'] = {group_id : samples for group_id, samples in config['groups'].items() if len(samples) == 1}
###

include: "genomes.smk"
# Define conveniance methods for accessing the fasta and gff files for downstream analysis.
# This allows you to provide your own reference or gff instead of automatically downloading from NCBI. 
def get_reference(wildcards):
    return rules.merge_annotations.output.fasta

def get_gff(wildcards):
    return rules.merge_annotations.output.gff

def get_contigs(wildcards):
    return rules.merge_annotations.output.contigs


include: "alignment.smk"
include: "postprocessing.smk"
include: "analysis.smk"


if config['_run_pipeline'] == 'build-ref':

    targets = [
        *rules.summarize_genomes.output,
        rules.make_index.output,
        rules.make_snpEff_db.output,
    ]

elif config['_run_pipeline'] == 'align':

    targets = list(expand(rules.markduplicates.output, sample = samples_list))
 
elif config['_run_pipeline'] == 'variants':
    
    targets = [
        rules.annotate_vcf.output,
        *expand(rules.summarize_abundances.output, sample = samples_list),
        *expand(rules.get_multimap_stats.output.stats, sample = samples_list)
    ]


# set the union of sample and group-level files as the target for the pipeline.
rule all:
    input: targets
        
