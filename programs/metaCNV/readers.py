import pyBigWig
import pyfaidx
import subprocess
import os
import numpy as np
import pandas as pd
import tempfile


def get_circular_ORI_distance(region, chrom_sizes, ORIs):
    
    region_pos = region.start + (region.end - region.start)/2
    chromlen = chrom_sizes[region.chr]
    ori_pos = ORIs[region.chr]

    direct_distance = region_pos - ori_pos
    if region_pos > ori_pos:
        # |~~~~~~~>O-----R~~~~~~~>|
        circular_distance = ori_pos + (chromlen - region_pos)
    else:
        # |<~~~~~~~R-----O<~~~~~~~|
        circular_distance = region_pos + (chromlen - ori_pos)

    # choose the minimum distance - 
    # direct or circular path between ORI and region elements
    # normalize so that the range is -0.5 to 0.5
    return min(abs(direct_distance), circular_distance)/chromlen*2 - 0.5



def get_coverage_from_bigwig(
    bigwig_path,
    region_gc,
):

    with pyBigWig.open(bigwig_path) as bw:
        
        region_coverages = np.zeros(len(region_gc))
        
        for i, (_, data) in enumerate(region_gc.iterrows()):
            
            chrom, start, end = data.chr, data.start, data.end
            try:
                region_coverages[i] = bw.stats(chrom, start, end)[0]
            except RuntimeError:
                region_coverages[i] = 0.
       
    return region_coverages



def get_gc_content_per_window(*,
    fasta_file, 
    outfile,
    window_size = 250,
):  
    
    index_file = fasta_file + '.fai'
    if not os.path.isfile(index_file):
        pyfaidx.index(fasta_file)

    gc_command = f"cut -f1,2 {index_file} | bedtools makewindows -g - -w {int(window_size)} | " \
                    f"bedtools nuc -bed - -fi {fasta_file} | cut -f1,2,3,5 | tail -n+2 > {outfile}"

    subprocess.run(
        gc_command, 
        shell = True, check = True
    )




def get_dataframe(*,
    fasta_file,
    bigwig_file,
    ori_file,
    window_size = 100,
    min_average_coverage = 5,
):
    
    with tempfile.NamedTemporaryFile() as windows_file:

        # 1. create windows and get GC content within each window
        get_gc_content_per_window(
                fasta_file = fasta_file,
                outfile = windows_file.name,
                window_size=window_size
            )
        
        regions = pd.read_csv(
                    windows_file.name, 
                    sep = '\t', 
                    names = ['chr', 'start', 'end', 'gc'], 
                    header = None
                )
        
        # 2. get coverage from bigwig file
        regions['coverage'] = get_coverage_from_bigwig(
            bigwig_path = bigwig_file,
            region_gc = regions,
        )  

    # 3. read contig sizes
    chrom_sizes = pd.read_csv(fasta_file + '.fai', 
                              sep = '\t', 
                              header = None,
                              names = ['chr','length'], 
                              usecols=[0,1], 
                              index_col = 0)\
                            ['length'].to_dict()

    # 4. read ORI positions
    oris = pd.read_csv(ori_file, header=None, sep = '\t',
                        names=['ORI','chr','start','end','cummulative_gc_skew','type'])                
    oris = oris[oris['type'] == 'origin']
    oris = (oris.set_index('ORI')['chr'] + window_size//2).to_dict()

    # 5. calculate distance to ORI for each region
    regions['ori_distance'] = regions.apply(
        lambda x : get_circular_ORI_distance(x, chrom_sizes, oris), 
        axis = 1
    )

    # 6. calculate the median abundance over each contig, this should be 
    #    approximately the expected coverage of non-CNV regions times 
    #    the abundance of the organism
    contig_abundances = regions[regions.coverage > 1]\
            .groupby('chr')['coverage'].median()

    # 7. filter out contigs with less than the minimum average coverage
    contig_abundances = contig_abundances[contig_abundances > min_average_coverage]

    regions = regions[ regions.chr.isin(contig_abundances.index) ]
    regions['exposure'] = regions.chr.map(contig_abundances)

    return regions