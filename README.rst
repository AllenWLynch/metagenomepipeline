SmilliePipeline - Variant calling from metagenomic samples
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Installation
------------

Create a new directory to house the pipeline and its associated files:

.. code-block:: bash

    $ mkdir pipeline && cd pipeline

And clone the Github repository:

.. code-block:: bash

    $ git clone https://github.com/AllenWLynch/metagenomepipeline.git

The pipeline only needs a few dependencies to run, all other dependencies
and tools are managed internally by snakemake. Create an environment for 
the pipeline using the "environment.yaml" file:

.. code-block:: bash

    $ conda env create -f environment.yaml
    $ conda activate smilliepipeline


Commands
--------

The pipeline comes with two command line interfactes, "smillie-pipeline",
and "pipeline-utils". The former is the main entrypoint for processing samples,
and the latter is used to manipulating configuration files. Make sure to 
make these programs executable, and add them to your path.

.. code-block:: bash

    $ smillie-pipeline --help
    
    usage: Command to run metagenomics pipeline [-h] {build-ref,align,count,variants} ...

    positional arguments:
        {build-ref,align,count,variants}
        build-ref           Download genomes and create indexes for various tools.
        align               Align reads to reference genome
        count               Count gene and genome coverages.
        variants            Call variants from aligned reads

The "smillie-pipeline" command has four subcommands, which launch nested pipelines.
For example, if you run "smillie-pipeline variants", this will also ensure that 
all the annotations are present by launching the "build-ref" pipeline, if needed. 
I suggest running these four commands in sequence to ensure that each previous
step performs as expected.

Each command takes the same arguments:

.. code-block:: bash

    $ ./smillie-pipeline variants --help
    
    usage: Command to run metagenomics pipeline variants [-h] --samples-config SAMPLES_CONFIG --genomes-config GENOMES_CONFIG --directory DIRECTORY
                                                        [--parameters-config PARAMETERS_CONFIG] [--cluster] [--cores CORES] [--snake-args ...]

    Pipeline arguments:
    --samples-config SAMPLES_CONFIG, -samples SAMPLES_CONFIG
                            Specifies the samples to process and their metadata.
    --genomes-config GENOMES_CONFIG, -genomes GENOMES_CONFIG
                            Configuration file for sequences and annotations
    --directory DIRECTORY, -dir DIRECTORY
                            Directory in which to save pipeline results. Providing a new directory starts a new pipeline run. If the directory already exists, the pipeline will
                            resume from the last completed step.
    --parameters-config PARAMETERS_CONFIG, -parameters PARAMETERS_CONFIG
                            Configuration file for parameters. Default: parameters.yaml

    Common arguments:
    --cluster             Run jobs on cluster instead of locally. This will use qsub to submit jobs to the scheduler.
    --cores CORES, -c CORES
                            Number of cores to use when running locally. Default: 1
    --snake-args ..., -s ...
                            Arguments to pass to snakemake. Every argument passed after the -s flag will be passed to snakemake. See the snakemake command line documentation
                            for more information.

**--samples-config**
Path to a yaml file which specifies the samples to process, and their metadata. That yaml file
should be in the format:

.. code-block:: yaml

    samples:
        sample_name_1:
            read1 : /path/to/R1.fastq.gz
            read2: /path/to/R2.fastq.gz
            is_paired: TRUE/FALSE
            is_trimmed: TRUE/FALSE
            metadata: 
                disease : IBD
                subject : S839
                cohort : PRISM
                key : value

        sample_name_2:
            ...

With an outer header "samples", under which each sample is listed. Each sample
has a unique name, with the following required fields:

- read1 : Path to the first read file. If the reads are single-end, then this is the only read file.
- read2 : Path to the second read file. If the reads are single-end, then this field is not required.
- is_paired : Whether the reads are paired-end or single-end. This is a boolean value, and can be TRUE/FALSE, True/False, true/false, 1/0, or yes/no.
- is_trimmed : Whether the reads have been trimmed. This is a boolean value, and can be TRUE/FALSE, True/False, true/false, 1/0, or yes/no.

The "metadata" field is optional, and can be used to specify any additional metadata about the sample.
You can provide any arbitrary key-value pairs.

**--genomes-config**
Path to a yaml file which specifies the genomes to use, and their
UHGG identifiers. This file should be in the format (for example):

.. code-block:: yaml

    genomes:
        p-dorei:
            GUT_ID: GUT_GENOME143505
            MGYG: MGYG-HGUT-02478
            species : Phocaeicola dorei
         f-prau:
            GUT_ID: GUT_GENOME140074
            MGYG: MGYG-HGUT-02272
            species : Faecalibacterium prausnitzii


With an outer header "genomes", under which each genome is listed. Each genome
has a unique name, and following required fields:

- GUT_ID : The UHGG identifier for the genome.
- MGYG : The MGYG identifier for the genome.
- species : The species name for the genome.

**--directory**
Specifies the directory in which to save the pipeline results.
If the directory already exists, the pipeline will resume from the last completed step for files within that directory, 
otherwise, the pipeline will start from the beginning.

**--parameters-config**
Specifies the path to a yaml file which contains the parameters
for the pipeline. By default, this is the "parameters.yaml" file in the pipeline directory,
which you can modify to change default parameters of the pipeline.

The next set of arguments control the pipeline runtime and adjust its behavior on the cluster.

**--cluster**
Tells the pipeline to run on the cluster instead of locally.
In this mode, each rule to be executed will be submitted as a qsub job.

**--cores**
Specifies the number of cores to use when running locally.

**--snake-args**
Every argument passed after the "--snake-args" flag will be passed to snakemake. See the snakemake command line documentation
for more information. Perhaps the most important argument is "--dryrun/-n", which will print out the commands that would be run
if the pipeline were executed to completion. I higly suggest running the pipeline with "--dryrun" first to ensure that the
pipeline is configured correctly before starting large jobs.

Another important argument is "--unlock", which will unlock the pipeline if it was previously interrupted. 

Finally, if you only want the pipeline to target a specific file/files, just append the filenames to the end of the command.

**Example**

Putting it all together, to test the full pipeline for some configuration on the Broad server, use:

.. code-block:: bash

    $ smillie-pipeline variants \
        -samples samples.yaml \
        -genomes genomes.yaml \
        --directory /path/to/results \
        --cluster \
        --snake-args --dryrun

Simply rerun this command without "--dryrun" to start processing samples.


Results
-------

The pipeline will save results in the directory with the following structure:

::

    ├── QC
    ├── analysis
    ├── benchmark
    ├── genomes
    ├── logs
    ├── processing


**analysis**: Where final, processed results are saved.
**analysis/all**: under "analysis/all", you will find summary statistics and results aggregated across samples. A completed pipeline run will produce:

- logratio-pileup.bw: the log10 fold ratio of unique primaries to multimapped reads across loci.
- multimap-pileup.bw: log10 of the number of multimapped reads across loci.
- primary-pileup.bw: log10 of the number of unique primaries across loci.
- multimap-stats.tsv: summary statistics for multimapped reads across loci, including the proportion of the genome which is covered with x fraction of multimapped reads.
- variants.bcf: final variant calls file

**analysis/samples**: results for each sample, with structure:

::

    ├── sample_name.bam
    ├── sample_name.bam.bai
    ├── coverage.bigwig
    ├── feature_counts
    │   ├── f-prau.tsv.summary
    │   └── p-dorei.tsv.summary
    ├── feature_counts.tsv
    └── multimapping
        ├── multimap-coverage.bigwig
        ├── multimap.bam
        ├── multimap.sorted-by-name.sam
        └── primary-coverage.bigwig

Where "feature_counts.tsv" contains the number of reads aligned to each element of each genome.

**genomes**: Stores annotation information, with a directory for each species' genome, as well as 
"all", which contains the concatenated annotations for all genomes. The information under "all" 
is used for alignment and variant calling, which works jointly across every species. The "all/contigs.tsv"
file summarizes the contigs associated with each genome and species:

::

    #genome #contig #species
    p-dorei GUT_GENOME143505_1      Phocaeicola-dorei
    f-prau  GUT_GENOME140074_1      Faecalibacterium-prausnitzii
    f-prau  GUT_GENOME140074_2      Faecalibacterium-prausnitzii
    f-prau  GUT_GENOME140074_3      Faecalibacterium-prausnitzii
    f-prau  GUT_GENOME140074_4      Faecalibacterium-prausnitzii

**benchmark**: Stores benchmarking information, including the time and memory usage of each rule.
**logs**: Stores the log files for each rule.
**QC** : Empty at the moment, but one could add FastQC or something if that was important.
**processing**: Stores intermediate files (SAM files, unfiltered VCFs, etc.), which are deleted when no longer needed.

Manipulating configuration files
--------------------------------

Often, it is easier to work with tabular files instead of yamls on the command line, so the "pipeline-utils" command 
can be used to convert between formats. The following line will convert a yaml file to a tabular format,
then back to yaml again.

.. code-block:: bash

    $ pipeline-utils to-df samples.yaml -s samples | pipeline-utils to-yaml -s samples - > samples.yaml

The "--subsection/-s" flag tells the command which top-level section from a yaml file to convert. In this case,
we open "samples.yaml", and convert the section with the heading "samples". When converting from tabular to yaml,
the "--subsection" flag is used to specify the top-level section to write to.


Pipeline architecture
---------------------

The pipeline is built using snakemake. Snakemake programs consist of a collection of "rules", 
which are instructions to turn some set of input files into some other set of output files by
running a command. Essentially, one provides snakemake with a list of files which you would 
like to have, then snakemake will figure out which rules need to be run to produce those files. 
Finally, snakemake determines the order to run those rules and handles job scheduling, resource
allocation, and dependency management.

The pipeline rules are laid out in the "rules" directory. Additionally, I have some
helper scripts in the "scripts" directory. The "query-gff" script is actually pretty useful
outside of the pipeline. When you invoke the command line, "smillie-pipeline" parses your
command line arguments and configurations, then starts a snakemake program. The program
first runs "rules/controller.smk", which sets the "target" of the pipeline - which files 
you would like to create. The pipeline subcommands (build-ref, variants, etc.) are actually
just aliases that instruct the controller to set the target to the appropriate files.