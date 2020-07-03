#!/usr/bin/env bash

## Download

mkdir -p GSE125384/fastq/Sample_EL4_rep1/
mkdir -p GSE125384/fastq/Sample_EL4_rep2/

wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR846/002/SRR8466982/SRR8466982.fastq.gz -O GSE125384/fastq/Sample_EL4_rep1/SRR8466982_R1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR846/003/SRR8466983/SRR8466983.fastq.gz -O GSE125384/fastq/Sample_EL4_rep2/SRR8466983_R1.fastq.gz

## Analysis

git clone git@gitlab.iric.ca:gendrop/dsptools.git

cp dsptools/rnaseq.GRCm38_M23.CC.conf rnaseq.conf
vi rnaseq.conf

./dsptools/run.sh -c rnaseq.conf -g qc-raw
./dsptools/run.sh -c rnaseq.conf -g star-rsem

./dsptools/run.sh -c rnaseq.conf rsem-annotate


## TPM filtering

(source ../filterRSEM.sh && filter genes && filter isoforms) | tee tpm.stats
