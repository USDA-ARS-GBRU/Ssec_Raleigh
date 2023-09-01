#!/bin/bash

if [[ $# -ne 4 ]]; then
    echo
    echo "Incorrect number of arguments"
    echo
    echo "Usage: "
    echo "   "$0" chrLengths chrNumToLength binSizeBp bcfgz"
    echo
    exit 1
fi


chrLengths=$1
chrNumToLength=$2
binSz=$3
bcfgz=$4

base=${bcfgz%.bcf.gz}
binSzKb=$(($binSz/1000))


awk -v binSz=$binSz '
BEGIN{
    FS=OFS="\t";
}
{   
    if (FNR==NR) {
        chrNum[$1] = FNR
    } else if ($1 !~ /^@/) {
        chr = chrNum[$1]
        print chr, $2, int($2/binSz)
    }
}
' $chrLengths <(bcftools view -H $bcfgz | cut -f1,2) > ${base}_chrPos${binSzKb}KBin.txt


awk -v binSz=$binSz '
BEGIN{
    FS=OFS="\t";
    maxChr=-999
}
{
    if(FNR==NR) {
        chr = $1
        if (chr > maxChr) {maxChr = chr}
        maxBin[chr] = int($2/binSz)
        for (i=0;i<=maxBin[chr];i++) {
            coverage[chr, i] = 0
        }
    } else {
        coverage[$1, $3]++
    }
}
END{
    for (chr=1; chr <= maxChr; chr++) {
        for (bin=0; bin <= maxBin[chr]; bin++) {
            print chr, bin, coverage[chr, bin]
        }
    }
}
' $chrNumToLength ${base}_chrPos${binSzKb}KBin.txt > ${base}_coverage${binSzKb}KBins.txt

