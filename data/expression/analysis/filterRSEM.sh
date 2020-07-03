#!/usr/bin/env bash

# Needs DSPTools: rsem-annotate to be executed to work

filter() {
    entry=$1  # genes or isoforms
    echo ${entry^^}
    echo

    s=`ls rsem | head -1`
    if [ -d rsem/$s/stranded ]
    then
        str=
    else
        str=un
    fi

    # Join TPM RSEM samples
    mkdir -p tpm
    gene_list=
    gene_list_prev=
    for file in `ls rsem/*/${str}stranded/*.$entry.results.cDNA.TPM.tsv`
    do
        echo $file
        gene_list=`cut -f 1 $file`
        if [ -z "$gene_list_prev" ]
        then
            cut -f 1,3 $file > tpm/$entry.tpm.tsv
        elif [ "$gene_list" == "$gene_list_prev" ]
        then
            cut -f 3 $file | paste -d "\t" tpm/$entry.tpm.tsv - > tpm/tempfile.tsv
            mv tpm/tempfile.tsv tpm/$entry.tpm.tsv
        else
            echo "ERROR" >&2
            break
        fi
        gene_list_prev=$gene_list
    done
    echo

    # Stats
    nentries=`cat rsem/*/${str}stranded/*.$entry.results.cDNA.TPM.tsv | grep -v transcript_id |\
        awk -F '\t' '$3>=0 {print $1}' | sort | uniq | wc -l`
    echo entries: $nentries

    nexpressed=`cat rsem/*/${str}stranded/*.$entry.results.cDNA.TPM.tsv | grep -v transcript_id |\
        awk -F '\t' '$3>0 {print $1}' | sort | uniq | wc -l`
    echo entries expressed: $nexpressed

    minval=`cat rsem/*/${str}stranded/*.$entry.results.cDNA.TPM.tsv | grep -v transcript_id |\
        awk -F '\t' '$3>0 {print $3}' | sort -n | head -1`
    echo entries minval: $minval
    echo

    # Filter
    awk -F '\t' '
            {sum=0;
             for (i=2; i <= NF; i++) sum += $i;
             sum /= (NF-1)}
            sum > 0 {t=$1; gsub("\\..*", "", t); print t"\t"sum}' \
        tpm/$entry.tpm.tsv > tpm/$entry.expressed.mean.tpm.txt

    nlines=`cat tpm/$entry.expressed.mean.tpm.txt | wc -l`
    echo mean expressed: $nlines

    minval=`cat tpm/$entry.expressed.mean.tpm.txt | cut -f 2 | sort -n |\
        awk -v n="$nlines" 'BEGIN {n=int(n/100)} NR >= n {print; exit}'`
    echo mean minval 99%: $minval

    cat tpm/$entry.expressed.mean.tpm.txt |\
        awk -v m="$minval" '$2>m' > tpm/expressed.99p.$entry.mean.tpm.txt
    nlines99=`cat tpm/expressed.99p.$entry.mean.tpm.txt | wc -l`
    echo mean expressed 99%: $nlines99
    echo

    awk -F '\t' '
            {n=NF-1;
             split("", a, ":");
             for (i=1; i <= n; i++) a[$i]=$(i+1);
             asort(a,b);
             (n % 2) ? med=b[int(n/2)+1] : med=(b[int(n/2)]+b[int(n/2)+1])/2}
            med > 0 {t=$1; gsub("\\..*", "", t); print t"\t"med}' \
        tpm/$entry.tpm.tsv > tpm/$entry.expressed.median.tpm.txt

    nlines=`cat tpm/$entry.expressed.median.tpm.txt | wc -l`
    echo median expressed: $nlines

    minval=`cat tpm/$entry.expressed.median.tpm.txt | cut -f 2 | sort -n |\
        awk -v n="$nlines" 'BEGIN {n=int(n/100)} NR >= n {print; exit}'`
    echo median minval 99%: $minval

    cat tpm/$entry.expressed.median.tpm.txt |\
        awk -v m="$minval" '$2>m' > tpm/expressed.99p.$entry.median.tpm.txt
    nlines99=`cat tpm/expressed.99p.$entry.median.tpm.txt | wc -l`
    echo median expressed 99%: $nlines99
    echo
}
