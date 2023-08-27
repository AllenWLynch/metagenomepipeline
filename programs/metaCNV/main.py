from readers import get_dataframe
from hmm import run_CMV_EM, regions_to_matrix
from sklearn.utils import check_random_state
from numpy import inf
import subprocess
from pandas import Series
import argparse
import logging

logger = logging.getLogger('metaCNV')
logger.setLevel(logging.INFO)


def get_ploidy(
    regions,
):

    model, scores = run_CMV_EM(regions)

    X, lengths = regions_to_matrix(regions)
    
    regions['ploidy'], regions['predicted_coverage'] = model.predict_coverage(X, lengths = lengths)

    return regions, model



def collapse_regions(
    regions
):
        
    chrom, start,end,ploidy = None,None,None,None

    for _, row in regions.iterrows():
            
        if chrom != row.chr or ploidy != row.ploidy or row.start != end:
            if chrom is not None:
                yield chrom, start, end, ploidy
            
            chrom, start,end,ploidy = row.chr, row.start, row.end, row.ploidy
        else:
            end = row.end
        
    yield chrom, start, end, ploidy



def save_results(
        regions, model,*,
        fasta_file,
        outprefix = 'metaCNV_analysis'
        ):
    
    logger.info(f'Writing results...')

    # 2. save the PTR data
    Series(model.get_PTR(), name = 'PTR').to_frame()\
        .reset_index().rename(columns = {'index' : 'genome'})\
        .to_csv(outprefix + '.PTR.tsv', index = False, sep = '\t')
    

    # 3. save a bedgraph matrix with all of the relevant information
    # produced by the analysis
    ploidy_matrix = outprefix + '.ploidy.bgmtx'
    regions.add_prefix('#', axis='columns')\
            .to_csv(ploidy_matrix, 
                    sep = '\t', 
                    index = False,
                    )
        
    subprocess.run(
        f'bgzip -f {ploidy_matrix} && tabix -p bed {ploidy_matrix}.gz', 
        shell= True, check=True
    )

    # save a bed file for use for variant calling with 
    # unannotated regions set to ploidy 1
    bedfile = outprefix + '.ploidy.bedgraph'
    index_file = fasta_file + '.fai'

    with open(bedfile, 'w') as f:
        for i, record in enumerate(collapse_regions(regions)):
                print(*record, sep = '\t', file = f, end = '\n')

    command = f"cut -f1,2 {index_file} | bedtools complement -i {bedfile} -g - | " \
        "awk -v OFS=\"\t\" '{print $0,\"1.0\"}' > tmp.bed && " \
        f"cat {bedfile} tmp.bed | sort -k1,1 -k2,2n | bgzip > {bedfile}.gz && " \
        f"tabix -p bed {bedfile}.gz && rm {bedfile}"
        

    subprocess.run(command, shell = True, check = True)
    


def main(*,
    fasta_file,
    bigwig_file,
    ori_file,
    outprefix = 'metaCNV_analysis',
    window_size = 100,
    min_average_coverage = 5,
    n_fits = 5
):

    save_results(
        *get_ploidy(
            get_dataframe(
                fasta_file = fasta_file,
                bigwig_file = bigwig_file,
                ori_file = ori_file,
                window_size = window_size,
                min_average_coverage = min_average_coverage,
            ),
        ),
        fasta_file = fasta_file,
        outprefix = outprefix,
        )


def get_parser():

    parser = argparse.ArgumentParser()
    parser.add_argument('--fasta','-fa', type = str, required = True)
    parser.add_argument('--bigwig','-bw', type = str, required = True)
    parser.add_argument('--ori','-ori', type = str, required = True)
    parser.add_argument('--outprefix','-o', type = str, required = True)
    parser.add_argument('--window_size','-ws', type = int, default = 100)
    
    return parser


if __name__ == '__main__':

    args = get_parser().parse_args()

    main(
        fasta_file = args.fasta,
        bigwig_file = args.bigwig,
        ori_file = args.ori,
        outprefix = args.outprefix,
        window_size = args.window_size,
    )