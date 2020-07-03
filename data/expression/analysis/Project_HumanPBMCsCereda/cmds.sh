#!/usr/bin/env bash

## Download

mkdir -p GSE106443/fastq/Sample_HealthyCtrl1/
mkdir -p GSE106443/fastq/Sample_HealthyCtrl2/
mkdir -p GSE106443/fastq/Sample_HealthyCtrl3/

wget -c ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR624/002/SRR6246152/SRR6246152_1.fastq.gz -O GSE106443/fastq/Sample_HealthyCtrl1/SRR6246152_R1.fastq.gz
wget -c ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR624/002/SRR6246152/SRR6246152_2.fastq.gz -O GSE106443/fastq/Sample_HealthyCtrl1/SRR6246152_R2.fastq.gz
wget -c ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR624/003/SRR6246153/SRR6246153_1.fastq.gz -O GSE106443/fastq/Sample_HealthyCtrl2/SRR6246153_R1.fastq.gz
wget -c ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR624/003/SRR6246153/SRR6246153_2.fastq.gz -O GSE106443/fastq/Sample_HealthyCtrl2/SRR6246153_R2.fastq.gz
wget -c ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR624/004/SRR6246154/SRR6246154_1.fastq.gz -O GSE106443/fastq/Sample_HealthyCtrl3/SRR6246154_R1.fastq.gz
wget -c ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR624/004/SRR6246154/SRR6246154_2.fastq.gz -O GSE106443/fastq/Sample_HealthyCtrl3/SRR6246154_R2.fastq.gz

## Analysis

git clone git@gitlab.iric.ca:gendrop/dsptools.git

cp dsptools/rnaseq.GRCh38_H32.CC.conf rnaseq.conf
vi rnaseq.conf

./dsptools/run.sh -c rnaseq.conf -g qc-raw
./dsptools/run.sh -c rnaseq.conf -g star-rsem

./dsptools/run.sh -c rnaseq.conf rsem-annotate

## TPM filtering

(source ../filterRSEM.sh && filter genes && filter isoforms) | tee tpm.stats
