
#configfile: "config.yaml"
"""
Config schema:

{
    # samples are individual sequencing experiments indexed by name
    "samples" : {
        "sample_name" : {
            "is_paired" : True,
            "read1" : "path/to/read1.fastq.gz",
            "read2" : "path/to/read2.fastq.gz", # if relevant
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
    
include: "alignment.smk"
include: "postprocessing.smk"
include: "analysis.smk"

# for each sample, collect these results
sample_level_files = [
    'htseq_counts/{sample}.htseq.counts',
    'coverage/{sample}.bamcoverage',
    'bams/{sample}.bam.bai',
]

sample_targets = [
        result
        for sample_result in sample_level_files
        for result in expand(sample_result, sample = config['samples'].keys())
    ]

# for each group collect these results
group_targets = expand(
    'vcfs/{group}.vcf.gz',
    group = config['groups'].keys()
)

# set the union of sample and group-level files as the target for the pipeline.
rule all:
    input: sample_targets + group_targets
        