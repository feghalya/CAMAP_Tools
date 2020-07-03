#!/usr/bin/env bash

## Mouse CT26 Perreault

p=MouseCT26_Perreault

xlsx2csv -d tab aau5516_Tables_S1_to_S18.xlsx -s 3 | tail -n +3 > aau5516_Table_S3.txt

cat aau5516_Table_S3.txt | tail -n +2 | wc -l
#^ 1875

# Get peptide list
cat aau5516_Table_S3.txt | tail -n +2 | cut -f 1 > ../../detected.peptides.${p}.txt

cd ../../..
