#!/usr/bin/env bash

## Download

mkdir -p GSE67174/fastq/Sample_BLCL_2.1/
mkdir -p GSE67174/fastq/Sample_BLCL_2.2/
mkdir -p GSE67174/fastq/Sample_BLCL_3/
mkdir -p GSE67174/fastq/Sample_BLCL_10/

wget -c ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR192/006/SRR1925276/SRR1925276_1.fastq.gz -O GSE67174/fastq/Sample_BLCL_2.1/SRR1925276_R1.fastq.gz
wget -c ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR192/006/SRR1925276/SRR1925276_2.fastq.gz -O GSE67174/fastq/Sample_BLCL_2.1/SRR1925276_R2.fastq.gz
wget -c ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR192/007/SRR1925277/SRR1925277_1.fastq.gz -O GSE67174/fastq/Sample_BLCL_2.1/SRR1925277_R1.fastq.gz
wget -c ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR192/007/SRR1925277/SRR1925277_2.fastq.gz -O GSE67174/fastq/Sample_BLCL_2.1/SRR1925277_R2.fastq.gz

wget -c ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR192/008/SRR1925278/SRR1925278_1.fastq.gz -O GSE67174/fastq/Sample_BLCL_2.2/SRR1925278_R1.fastq.gz
wget -c ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR192/008/SRR1925278/SRR1925278_2.fastq.gz -O GSE67174/fastq/Sample_BLCL_2.2/SRR1925278_R2.fastq.gz
wget -c ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR192/009/SRR1925279/SRR1925279_1.fastq.gz -O GSE67174/fastq/Sample_BLCL_2.2/SRR1925279_R1.fastq.gz
wget -c ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR192/009/SRR1925279/SRR1925279_2.fastq.gz -O GSE67174/fastq/Sample_BLCL_2.2/SRR1925279_R2.fastq.gz

wget -c ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR192/000/SRR1925280/SRR1925280_1.fastq.gz -O GSE67174/fastq/Sample_BLCL_3/SRR1925280_R1.fastq.gz
wget -c ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR192/000/SRR1925280/SRR1925280_2.fastq.gz -O GSE67174/fastq/Sample_BLCL_3/SRR1925280_R2.fastq.gz
wget -c ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR192/001/SRR1925281/SRR1925281_1.fastq.gz -O GSE67174/fastq/Sample_BLCL_3/SRR1925281_R1.fastq.gz
wget -c ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR192/001/SRR1925281/SRR1925281_2.fastq.gz -O GSE67174/fastq/Sample_BLCL_3/SRR1925281_R2.fastq.gz

wget -c ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR192/002/SRR1925282/SRR1925282_1.fastq.gz -O GSE67174/fastq/Sample_BLCL_10/SRR1925282_R1.fastq.gz
wget -c ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR192/002/SRR1925282/SRR1925282_2.fastq.gz -O GSE67174/fastq/Sample_BLCL_10/SRR1925282_R2.fastq.gz
wget -c ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR192/003/SRR1925283/SRR1925283_1.fastq.gz -O GSE67174/fastq/Sample_BLCL_10/SRR1925283_R1.fastq.gz
wget -c ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR192/003/SRR1925283/SRR1925283_2.fastq.gz -O GSE67174/fastq/Sample_BLCL_10/SRR1925283_R2.fastq.gz


## Analysis

git clone git@gitlab.iric.ca:gendrop/dsptools.git

cp rnaseq.GRCh38_H32.CC.conf rnaseq.conf
vi rnaseq.conf

./dsptools/run.sh -c rnaseq.conf -g qc-raw
./dsptools/run.sh -c rnaseq.conf -g star-rsem

./dsptools/run.sh -c rnaseq.conf rsem-annotate

## TPM filtering

(source ../filterRSEM.sh && filter genes && filter isoforms) | tee tpm.stats
