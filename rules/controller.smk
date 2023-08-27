
genomes_list = [metadata['reference'] for species, metadata in config['genomes'].items()]
samples_list = list(config['samples'].keys())

def double_on_failure(base_resources):
    def _double_on_failure(wildcards, attempt):
        return base_resources * attempt

    return _double_on_failure

include: "genomes.smk"
# Define conveniance methods for accessing the fasta and gff files for downstream analysis.
# This allows you to provide your own reference or gff instead of automatically downloading from NCBI. 
def get_reference(wildcards):
    return rules.merge_annotations.output.fasta

def get_gff(wildcards):
    return rules.merge_annotations.output.gff


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
    ]


# set the union of sample and group-level files as the target for the pipeline.
rule all:
    input: targets
        
