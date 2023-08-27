import argparse
import sys
import faidx

def get_most_common_allele(ref, alts, allele_counts, total_count):

    allele_counts = dict(zip( alts.split(','), map(int, allele_counts.split(',')) ))
    allele_counts[ref] = total_count - sum(allele_counts.values())

    return max(allele_counts)


if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument('input', type = str)
    parser.add_argument('output', type = argparse.FileType('w'), default=sys.stdout)
    parser.add_argument('--reference','-ref', type = str)
    parser.add_argument('--region','-r', type = str)
    args = parser.parse_args()

    chrom, start, end = args.region.split('[:|-]')

    with faidx.Fasta(args.reference) as fa:
        reference_seq = fa[chrom][start:end].seq

    for chrom, pos, *allele_counts in args.input.readlines():
        
        consensus_allele = get_most_common_allele(*allele_counts)

        reference_seq[pos - start:pos-start+len(consensus_allele)] = consensus_allele

    
    print(reference_seq, file = args.output)
    

    