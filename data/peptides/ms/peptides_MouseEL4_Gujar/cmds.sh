#!/usr/bin/env bash

## Mouse EL4 Gujar

p=MouseEL4_Gujar

cat pr6b00971_si_001_EL4.txt |\
    awk -F '\t' '$1 == "Search" {next} 1' | cut -f 4 | sort | uniq | wc -l
#^ 3757

# Get peptide list
cat pr6b00971_si_001_EL4.txt |\
    awk -F '\t' '$1 == "Search" {next} 1' | cut -f 4 | sort | uniq > ../../detected.peptides.${p}.txt
