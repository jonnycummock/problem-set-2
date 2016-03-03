#! /usr/bin/env bash


# Question 1 Use BEDtools intersect to identify the size of the largest overlap between CTCF and H3K4me3 locations.

datasets=/Users/jonny/data-sets
tfbs_bed=$datasets/bed/encode.tfbs.chr22.bed.gz
histone_bed=$datasets/bed/encode.h3k4me3.hela.chr22.bed.gz

ansQ1=$(gzcat $tfbs_bed | awk '$4 == "CTCF"' > ctcf-peaks.bed

bedtools intersect -a ctcf-peaks.bed -b $histone_bed -wo \
    | awk '{print $NF}' \
    | sort -nr | head -n1)

echo "answer-1:" $ansQ1

# -wo means to add show A and B and the in C the overlap
# $7 refers to 7th column. $NF refers to the last column, wherever it is


# Question 2 Use BEDtools to calculate the GC content of nucleotides 19,000,000 to 19,000,500 on chr22 of hg19 genome build. Report the GC content as a fraction (e.g., 0.50).

fasta=$datasets/fasta/hg19.chr22.fa

interval_bed=$datasets/bed/test.bed

ansQ2=$(echo -e "chr22\t19000000\t19000500" > $datasets/bed/test.bed

bedtools nuc -fi $fasta -bed $interval_bed \
    | cut -f5 \
    | grep -v "pct")

echo "answer-2:" $ansQ2


# you have to build a bedfile that has coordinates and feed it into nuc
# to build the bed file, you are going to "echo -e" it

# the output of the map is mean signal 

# Question 3 Use BEDtools to identify the length of the CTCF ChIP-seq peak (i.e., interval) that has the largest mean signal in ctcf.hela.chr22.bg.gz.

# the output of the map tool is mean signal

ansQ3=$(bedtools map -a $datasets/bed/ctcf-peaks.bed -b $datasets/bedtools/ctcf.hela.chr22.bg -c 4 -o mean \
    | sort -k5nr \
    | head -n1 \
    | awk '{print $3 - $2 }')

echo "answer-3:" $ansQ3

# Question 4 Use BEDtools to identify the gene promoter (defined as 1000 bp upstream of a TSS) 
# with the highest median signal in ctcf.hela.chr22.bg.gz. Report the gene name (e.g., 'ABC123')

signal=$datasets/bedtools/ctcf.hela.chr22.bg
tss=$datasets/bed/tss.hg19.chr22.bed
genome=$datasets/genome/hg19.genome

bedtools flank -i $tss -g $genome -l 1000 -r 0 -s \
    | bedtools  sort -i - >  promoters.bed

ansQ4=$(bedtools map -a promoters.bed -b $signal -c 4 -o median \
    | sort -k7nr \
    | head -n1 \
    | cut -f4)

echo "answer-4:" $ansQ4

# given a tss file, how do you find promoters a thousand bases upstream in
# tss? with bedtools. you might be able to use slop, however bedtools
# flank is going to be your best bet. look up bedtools flank. there is
# also a -s to pay attention to strand (on a + strand the promoter is on


# the left side). so now you have promoters, intervals, and now you can
# compare that using the map tool. the summary statistic is median. for
# the highest median signal, give him the gene name

# Question 5 Use BEDtools to identify the longest interval on chr22 that is not 
# covered by genes.hg19.bed.gz. Report the interval like chr1:100-500.

genes=$datasets/bed/genes.hg19.bed

ansQ5=$(bedtools complement -i $genes -g $genome \
    | bedtools sort -i - \
    | grep "chr22" \
    | awk '{print $1, $2, $3, $3 - $2}' \
    | sort -k4nr \
    | head -n1 \
    | awk 'BEGIN {OFS=""} {print $1, ":", $2, "-", $3}')

echo "answer-5:" $ansQ5

# the tool you are going to want to use for this is bedtools complement.
# given input intervals, it will give you back all the intervals not
# covered by the input. you have to give it a genome file. you will get a
# bunch of new intervals and you then have to find the biggest one

# Question 6 (extra credit) Use one or more BEDtools that we haven't covered in class. Be creative.

tss=$datasets/bedtools/class5/tss.bed
tss2=$datasets/bedtools/class5/tss2.bed

sortBed -i $tss > tss.sorted.bed
sortBed -i $tss2 > tss2.sorted.bed

ansQ6=$(bedtools jaccard -a tss.sorted.bed -b tss2.sorted.bed)

echo "answer-6:" $ansQ6
# look up the full list of bedtools sub-commands to get ideas
