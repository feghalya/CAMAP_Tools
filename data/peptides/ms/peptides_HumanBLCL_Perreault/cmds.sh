#!/usr/bin/env bash

## Human B-LCL Perreault

p=HumanBLCL_Perreault

# Get peptide list
cat raw_files_corrected/*.tsv |\
    awk -F '\t' '(length($6) >= 8) && (length($6) <= 11) && ($(NF-1) <= 1) {print $6}' |\
    sort | uniq > ../../detected.peptides.${p}.txt

cat JCI88590.sdt2.tsv | tail -n +2 | cut -f 1 |\
    sort | uniq > ../../peptides.paper.${p}.txt

cat BLCLdataset.csv | tail -n +2 | awk -F ',' '$2 == "SOURCE" {print $5}' |\
    sort | uniq > ../../peptides.all.${p}.txt
cat BLCLdataset.csv | tail -n +2 | awk -F ',' '$3 == "test" && $2 == "SOURCE" {print $5}' |\
    sort | uniq > ../../peptides.test.${p}.txt
cat BLCLdataset.csv | tail -n +2 | awk -F ',' '$3 == "validation" && $2 == "SOURCE" {print $5}' |\
    sort | uniq > ../../peptides.validation.${p}.txt
cat BLCLdataset.csv | tail -n +2 | awk -F ',' '$3 == "train" && $2 == "SOURCE" {print $5}' |\
    sort | uniq > ../../peptides.train.${p}.txt

wc -l ../../*${p}*
#^
#  50615 detected.peptides.txt
#  25271 detected.peptides.paper.txt
#   3574 detected.peptides.test.txt
#   3603 detected.peptides.validation.txt
#   9548 detected.peptides.train.txt
#  14660 detected.peptides.all.txt

# Verify that selected peptides for training are containg in full dataset
awk -F '\t' 'NR==FNR {a[$1]++; next} ! a[$1]' \
    ../../detected.peptides.${p}.txt \
    ../../peptides.all.${p}.txt | wc -l
#^ 0
# note that when lowering FDR cutoff ($(NF-1) <= 1) we start having mismatches
