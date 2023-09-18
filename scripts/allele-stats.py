#!/usr/bin/env python3

import sys
import argparse
from math import log
from collections import Counter

def alleles_equal(a1, a2):

    assert len(a1) == len(a2)

    for n1, n2 in zip(a1, a2):
        if not n1 == 'N' and not n2 == 'N':
            if not n1 == n2:
                return False

    return True


def get_num_alleles(sample_alleles):

    allele_list = []
    n_samples = 0
    
    for sample_allele in sample_alleles:
        
        if not sample_allele.count('N') == len(sample_allele):
            n_samples+=1
            
            if len(allele_list) == 0:
                allele_list.append(sample_allele)
            else:
                for alt_allele in allele_list:
                    if alleles_equal(sample_allele, alt_allele):
                        break
                else:
                    allele_list.append(sample_allele)

    return len(allele_list)/log(n_samples), allele_list, n_samples



def get_consensus(unique_alleles):
    consensus_allele = ''.join(
        [Counter(pos).most_common(1)[0][0] 
         for pos in list(zip(*unique_alleles))
        ])

    return consensus_allele



def get_num_variants_distribution(unique_alleles, consensus_allele):

    diffcounts = []
    for allele in unique_alleles:
        diffcounts.append(
            sum( [a1 != a2 for a1,a2 in zip(allele, consensus_allele) if not a1=='N' and not a2=='N' ])
        )
    
    return diffcounts


def main(allele_list, output):

    theta_hat, alleles, n_samples = get_num_alleles(
                map(lambda x : x.strip(), allele_list.readlines())
            )

    consensus_allele = get_consensus(alleles)

    diffcounts = get_num_variants_distribution(alleles, consensus_allele)

    print(
        theta_hat, len(alleles), n_samples, consensus_allele, ','.join(map(str, diffcounts)),
        sep = '\t', end = '\n', file = output
    )
    

if __name__=="__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument('allele_list', type = argparse.FileType('r'), default = sys.stdin)
    parser.add_argument('--output', type = argparse.FileType('w'), default = sys.stdout)
    args = parser.parse_args()

    main(args.allele_list, args.output)
