#!/usr/bin/env bash

## Human B721.221 Wu

p=HumanB721_Wu

cat hits_16_* | grep -v '^allele' | cut -d ' ' -f 3 | sort | uniq | wc -l
#^ 21499

# Get peptide list
cat hits_16_* | grep -v '^allele' | cut -d ' ' -f 3 | sort | uniq > ../../detected.peptides.${p}.txt
