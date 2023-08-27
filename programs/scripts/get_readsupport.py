import argparse
import sys

def get_allele_counts(chrom, pos, ref, alts, af, total_count):

    alts = alts.split(',')
    assert all([len(alt) == 1 for alt in alts])
    assert len(ref) == 1

    allele_counts = dict(zip( alts, map(int, allele_counts.split(',')) ))
    allele_counts[ref] = total_count - sum(allele_counts.values())

    return [chrom, pos, *list(allele_counts.values())]


if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument('input', type = argparse.FileType('r'), default=sys.stdin)
    parser.add_argument('output', type = argparse.FileType('w'), default=sys.stdout)

    args = parser.parse_args()

    #add header
    print('#' + '\t#'.join(['chrom', 'pos', 'A', 'T', 'G', 'C', '.']), file = args.output)

    for line in args.input.readlines():
        
        if not line.startswith('#'):

            record = line.strip().split('\t')
            print(
                *get_allele_counts(*record),
                sep = '\t', 
                file = args.output
            )