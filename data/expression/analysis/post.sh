#!/usr/bin/env bash

for d in Project_*
do
    pushd $d
    d=${d/Project_/}
    for f in tpm/expressed.99p*
    do
        cp $f ../../`basename ${f/expressed.99p./} | sed 's/.txt$//'`.99p.$d.tsv
    done
    popd
done
