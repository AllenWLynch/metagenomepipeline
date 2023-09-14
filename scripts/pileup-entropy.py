#!/usr/bin/env python3

import sys
import argparse
from collections import Counter
from math import log2

def main(pileup, model):

    for line in pileup:
        
        fields = line.strip().split('\t')

        ref_base = fields[2].upper()
        check_groups = list('ATGC')

        basecalls = Counter( fields[4].upper().replace(',',ref_base).replace('.', ref_base) )
        num_reads = sum([basecalls[v] for v in check_groups])

        if model == 'bases':
            entropy = -sum(
                [basecalls[base]/num_reads * log2(basecalls[base]/num_reads) for base in check_groups if basecalls[base] > 0]
            )
        elif model == 'major-minor':
            
            numref = basecalls[ref_base]
            entropy = \
                -(numref/num_reads * log2(numref/num_reads) if numref > 0 else 0) \
                -( (1 - numref/num_reads) * log2(1 - numref/num_reads) if numref < num_reads else 0)
            
        else:
            raise ValueError()


        basecall_summary = ' '.join([f'{(base if base == ref_base else base.lower())}:{basecalls[base]}' for base in check_groups])

        print(fields[0], int(fields[1])-1, fields[1], basecall_summary, entropy, 
              sep = '\t', end = '\n', file = sys.stdout)

if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('pileup', type = argparse.FileType('r'), default = sys.stdin)
    parser.add_argument('--model', choices=['major-minor','bases'], default = 'major-minor')
    args = parser.parse_args()

    main(args.pileup, args.model)


        
