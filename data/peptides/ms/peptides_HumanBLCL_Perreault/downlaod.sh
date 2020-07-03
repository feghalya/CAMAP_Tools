#!/usr/bin/env bash

# paper dataset containing HLA-A and HLA-B peptides
wget https://dm5migu4zj3pb.cloudfront.net/manuscripts/88000/88590/JCI88590.sdt1-6.xlsx
xlsx2csv -d tab JCI88590.sdt1-6.xlsx -s 2 | tail -n +2 > JCI88590.sdt2.tsv

# same MS dataset containing HLA-A, HLA-B and HLA-C peptides
mkdir raw_files
cd raw_files
rsync -au feghalya@cluster.iric.ca:/u/daoudat/py/imm_prediction/ImmPred/ImmPred/data/raw_files/Pat* .
cd ..
mkdir raw_files_corrected
for f in raw_files/*
do
    f=`basename $f`
    cat raw_files/$f | tr -d '\r' |\
        awk '(NR == 1) && (substr($0, 1, 7) == "Peptide") {repair=1; print ","$0; next}
             repair {print "NA,"$0} ! repair {print}' |\
        awk -F ',' '{split($0, a, "\"");
                     for (i=1; i<=length(a); i++) {
		         if (i%2) gsub(",", "\t", a[i]);
		         printf a[i]"\""
                     };
                     printf "\n"
                    }' |\
        sed 's/"$//' > raw_files_corrected/${f/.csv/.tsv}
done
