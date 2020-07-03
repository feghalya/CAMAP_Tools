#!/usr/bin/env bash

## Human PBMCs Gujar

p=HumanPBMCs_Gujar

cat pr6b00971_si_002_PBMC* |\
    awk -F '\t' '$1 == "Search" {next} 1' | cut -f 2 | sort | uniq | wc -l
#^ 10567

cat pr6b00971_si_002_PBMC* |\
    awk -F '\t' '$1 == "Search" {next} $15 < 1250 || $17 < 1250 || $19 < 1250 || $21 < 1250' |\
    cut -f 2 | sort | uniq | wc -l
#^ 7422

cat pr6b00971_si_002_PBMC* |\
    awk -F '\t' '$1 == "Search" {next} $15 < 50 || $17 < 50 || $19 < 50 || $21 < 50' |\
    cut -f 2 | sort | uniq | wc -l
#^ 3644

cat pr6b00971_si_002_PBMC* |\
    awk -F '\t' '$1 == "Search" {next} $16 < 1 || $18 < 1 || $20 < 1 || $22 < 1' |\
    cut -f 2 | sort | uniq | wc -l
#^ 6556

# Get peptide list

# Do not filter on NetMHC prediction for now
cat pr6b00971_si_002_PBMC* |\
    awk -F '\t' '$1 == "Search" {next} 1' | cut -f 2 | sort | uniq > ../../detected.peptides.${p}.txt
# Filter based on rank as recommended by NetMHC
cat pr6b00971_si_002_PBMC* |\
    awk -F '\t' '$1 == "Search" {next} $16 < 1 || $18 < 1 || $20 < 1 || $22 < 1' |\
    cut -f 2 | sort | uniq > ../../peptides.rank1.${p}.txt
# Filter based on binding affinity
cat pr6b00971_si_002_PBMC* |\
    awk -F '\t' '$1 == "Search" {next} $15 < 50 || $17 < 50 || $19 < 50 || $21 < 50' |\
    cut -f 2 | sort | uniq > ../../peptides.50nM.${p}.txt
