#!/bin/bash

fasta=$1
chromsizes=$2
outprefix=$3
windowsize=10000

if [ -z "$chromsizes" ] || [ -z "$outprefix" ] || [ -z "$fasta" ]; then
    echo "Usage: $0 fasta chromsizes outprefix"
    exit 1
fi

intermediate_output=$outprefix

bedtools makewindows -g $chromsizes -w $windowsize -s 100 | \
    bedtools nuc -bed - -fi $fasta | tail -n+2 | \
    awk -v OFS="\t" -v windowsize=$windowsize \
        '{skew=($8-$7)/($8+$7); cumsum[$1]+=skew; print $1,$2+windowsize/2-50,$2+windowsize/2+50,skew,cumsum[$1]}' | \
    bgzip > $intermediate_output && \
    tabix -s1 -b2 -e3 -0 $intermediate_output

gunzip -c $intermediate_output | awk -v OFS="\t" '(!max[$1] || $5>=max[$1]){max[$1]=$5; maxrow[$1]=$0} (!min[$1] || $5<=min[$1]){min[$1]=$5; minrow[$1]=$0} END{for (chrom in maxrow) print maxrow[chrom],"terminus"; for (chrom in minrow) print minrow[chrom],"origin"}' | \
    sort -k1,1 -k2,2n #> $outprefix.origin_terminus.bed