#!/usr/bin/env bash

## Download

mkdir -p GSE111092/fastq/Sample_CT26/

wget -c ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR678/005/SRR6781405/SRR6781405_1.fastq.gz -O GSE111092/fastq/Sample_CT26/SRR6781405_R1.fastq.gz
wget -c ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR678/005/SRR6781405/SRR6781405_2.fastq.gz -O GSE111092/fastq/Sample_CT26/SRR6781405_R2.fastq.gz

## Analysis

git clone git@gitlab.iric.ca:gendrop/dsptools.git

cp dsptools/rnaseq.GRCm38_M23.CC.conf rnaseq.conf
vi rnaseq.conf

./dsptools/run.sh -c rnaseq.conf -g qc-raw
./dsptools/run.sh -c rnaseq.conf -g star-rsem

./dsptools/run.sh -c rnaseq.conf rsem-annotate

## TPM filtering

(source ../filterRSEM.sh && filter genes && filter isoforms) | tee tpm.stats
