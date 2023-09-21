import pandas as pd
from collections import Counter
import argparse
import sys

def main(*,contigs_file, chromsizes_file, bedgraph, outfile):

    contigs = pd.read_csv(contigs_file, sep = '\t')
    chromsizes = pd.read_csv(chromsizes_file, sep = '\t', header = None, names = ['#contig','#size'])

    genome_sizes = contigs.merge( chromsizes, on = '#contig' ).groupby('#genome')['#size'].sum().to_dict()
    
    contig_genome_map = dict(zip( contigs['#contig'].values, contigs['#genome'].values ))    

    thresholds = [1,3,10,100,1000,10000,100000]
    len_trackers = { genome : Counter() for genome in contig_genome_map.values() }

    for line in bedgraph:
        if not line.startswith('#'):
            contig,start,end,count = line.strip().split('\t')
            start,end,count = int(start),int(end),int(count)
            
            for threshold in thresholds:
                if count >= threshold:
                    len_trackers[ contig_genome_map[contig] ][threshold] += end - start

    

    print('#genome','#genome_size',*[f'#gt_{threshold}x_multimap_coverage' for threshold in thresholds], sep = '\t', file = outfile)

    for genome, tresholds in len_trackers.items():

        print(
            genome, genome_sizes[genome], *[int(tresholds[threshold]) for threshold in thresholds], sep = '\t', file = outfile,
        )

if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('--contigs-file', required = True)
    parser.add_argument('--chromsizes-file', required = True)
    parser.add_argument('--bedgraph-file','-bg', default = sys.stdin, type=argparse.FileType('r'))
    parser.add_argument('--output-file', '-o', default=sys.stdout, type=argparse.FileType('w'))

    args = parser.parse_args()

    main(
        contigs_file = args.contigs_file,
        chromsizes_file = args.chromsizes_file,
        bedgraph = args.bedgraph_file,
        outfile = args.output_file,
    )