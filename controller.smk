"""
Config schema:

{
    # samples are individual sequencing experiments indexed by name
    "samples" : {
        "sample_name" : {
            "is_paired" : True,
            "read1" : "path/to/read1.fastq.gz",
            "read2" : "path/to/read2.fastq.gz", # if relevant
            "metadata" :
                "arbitrary_data" : data,
        },
    },

    # groups are groups of samples belonging to the same subject 
    # - variants will be called toegher for grouped samples
    "groups" : {
        "group_name" : [
            "sample_name1", "sample_name2"
        ],
    },

    "index" : "path/to/bt2_index",
}
"""


def double_on_failure(base_resources):
    def _double_on_failure(wildcards, attempt):
        return base_resources * attempt

    return _double_on_failure


genomes_list = [metadata['reference'] for species, metadata in config['genomes'].items()]
samples_list = list(config['samples'].keys())

include: "genomes.smk"
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
        
