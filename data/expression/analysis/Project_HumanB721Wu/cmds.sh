#!/usr/bin/env bash

## Download

mkdir -p GSE93315/fastq/Sample_B721.221_A2902/
mkdir -p GSE93315/fastq/Sample_B721.221_A5101/
mkdir -p GSE93315/fastq/Sample_B721.221_A5401/
mkdir -p GSE93315/fastq/Sample_B721.221_A5701/

wget -c ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR516/007/SRR5163127/SRR5163127_1.fastq.gz -O GSE93315/fastq/Sample_B721.221_A2902/SRR5163127_R1.fastq.gz
wget -c ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR516/007/SRR5163127/SRR5163127_2.fastq.gz -O GSE93315/fastq/Sample_B721.221_A2902/SRR5163127_R2.fastq.gz
wget -c ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR516/008/SRR5163128/SRR5163128_1.fastq.gz -O GSE93315/fastq/Sample_B721.221_A5101/SRR5163128_R1.fastq.gz
wget -c ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR516/008/SRR5163128/SRR5163128_2.fastq.gz -O GSE93315/fastq/Sample_B721.221_A5101/SRR5163128_R2.fastq.gz
wget -c ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR516/009/SRR5163129/SRR5163129_1.fastq.gz -O GSE93315/fastq/Sample_B721.221_A5401/SRR5163129_R1.fastq.gz
wget -c ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR516/009/SRR5163129/SRR5163129_2.fastq.gz -O GSE93315/fastq/Sample_B721.221_A5401/SRR5163129_R2.fastq.gz
wget -c ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR516/000/SRR5163130/SRR5163130_1.fastq.gz -O GSE93315/fastq/Sample_B721.221_A5701/SRR5163130_R1.fastq.gz
wget -c ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR516/000/SRR5163130/SRR5163130_2.fastq.gz -O GSE93315/fastq/Sample_B721.221_A5701/SRR5163130_R2.fastq.gz


## Analysis

git clone git@gitlab.iric.ca:gendrop/dsptools.git

cp dsptools/rnaseq.GRCh38_H32.CC.conf rnaseq.conf
vi rnaseq.conf

./dsptools/run.sh -c rnaseq.conf -g qc-raw
./dsptools/run.sh -c rnaseq.conf -g star-rsem

./dsptools/run.sh -c rnaseq.conf rsem-annotate

## TPM filtering

(source ../filterRSEM.sh && filter genes && filter isoforms) | tee tpm.stats
