import argparse
import sys
import faidx
import tabix

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

    with tabix.open(args.input.name) as tab:

        for _, pos, *allele_counts in map(
            lambda x : x.strip().split('\t'), 
            tab.query(chrom, start, end)
        ):

            def max_tie_reference(allele_counts, refnuc):
                pass

            new_allele = 'ATGC-'[allele_counts.index( 
                max_tie_reference(allele_counts, reference_seq[pos - start]) 
            )]
        
            reference_seq[pos - start] = new_allele

    print(reference_seq)

    

    