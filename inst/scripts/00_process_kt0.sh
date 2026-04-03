#!/bin/bash

# script to be run out of extdata folder
cd kt/inst/extdata
mkdir -p rename

# For each psi-blastp / blastp search
# the complete sequence output for each seed
# (filename) is in the hits/ folder
# rename the files to be one line for CSV import
for file in $(ls nr_hits/)
do
	echo $file
	seqkit replace -w 0 -p "^" -r "#$file@" -s hits/$file \
	| tr -sd '\n' - \
	| sed 's/#/\t/g' - \
	| sed 's/@/\t/g' - \
	| sed 's/>/\n>/g' - > rename/$file

done

cat hits/*   > kt0.raw.fa
cat rename/* > kt0.raw.tsv

# These are manually imported and renamed in
# a spreadsheet to yield kt0.preclust.fa
usearch --cluster_fast kt0.preclust.fa \
  -id 0.90 \
  -centroids ../kt1/kt1.fa \
  -uc ../kt1/kt1.id90.uc

# All vs. All Pairwise Alignment
#query target id alnlen mism opens qlo qhi tlo thi evalue bits
usearch -allpairs_local ../kt1/kt1.fa \
        -id 0.25 -evalue 0.1 \
        -blast6out ../kt1/kt1.aln