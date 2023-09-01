#!/bin/bash

if [[ $# -ne 5 ]]; then
    echo
    echo "Incorrect number of arguments"
    echo
    echo "Usage: "
    echo "   "$0" species binSizeBp minMAPQ bamDir bam"
    echo
    exit 1
fi


species=$1
binSz=$2
minQ=$3
bamDir=$4
bam=$5

base=${bam%.bam}
binSzKb=$(($binSz/1000))

# get the chromosome lengths (if you havent already)
[ -f ${species}.chrLengths.txt ] || samtools view -H ${bamDir}/${bam} | grep "^@SQ" | cut -f 2,3 | grep "^SN:" | tr "\t" ":" | cut -d : -f 2,4 | tr ":" "\t" | sort -V -k 1,1 > ${species}.chrLengths.txt


[ -f ${species}.chrNumToLength.txt ] || awk '
BEGIN{
    FS=OFS="\t";
}
{
    print NR,$2
}
' ${species}.chrLengths.txt > ${species}.chrNumToLength.txt


[ -f ${species}.chrNumToChrName.txt ] || awk '
BEGIN{
    FS=OFS="\t";
}
{
    print NR,$1
}
' ${species}.chrLengths.txt > ${species}.chrNumToChrName.txt


awk -v binSz=$binSz '
BEGIN{
    FS=OFS="\t";
}
{   
    if (FNR==NR) {
        chrNum[$1] = FNR
    } else if ($1 !~ /^@/) {
        chr = chrNum[$3]
        print chr, $4, int($4/binSz)
    }
}
' ${species}.chrLengths.txt <(samtools view -f PROPER_PAIR,READ1 -F UNMAP,SECONDARY,DUP -q $minQ ${bamDir}/${bam}) > ${base}_chrPos${binSzKb}KBin.txt


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
' ${species}.chrNumToLength.txt ${base}_chrPos${binSzKb}KBin.txt > ${base}_coverage${binSzKb}KBins.txt

