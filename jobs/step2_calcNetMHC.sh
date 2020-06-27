#!/usr/bin/env bash

name=calc_netMHC
date=`date +%Yy%mm%dd_%Hh%Mm%Ss`

njobs=0

function runPBS {
    genome=$1
    fullname=${name}_${genome}
    echo $fullname

    netmhc_dir=NetMHCpan-4.0a

    pushd output/allPeptides_$genome

    if [ ! -d $netmhc_dir/peptides ]
    then
        mkdir $netmhc_dir/peptides
        # remove selenocysteines not understood by netMHCpan
        time cat proteome.wide.peptides.txt | grep -v U |\
            awk -v dir=$netmhc_dir/peptides '{
                    letter=substr($0, 1, 1);
                    len=length($0);
                    print $0 >> dir"/"letter"."len}'
    fi

    if [[ $genome == "GRCh"* ]]
    then
        alleles=`cat ../../data/alleles/Human.tsv | cut -f 1 | tail -n +2 | head -n -1`
    else
        alleles=`cat ../../data/alleles/Mouse.tsv | cut -f 1 | tail -n +2`
    fi

    mkdir -p  $netmhc_dir/jobs
    c=0
    for f in `ls $netmhc_dir/peptides`
    do
        allele_count=0
        cur_alleles=
        for all in $alleles
        do
            if [[ $all == "HLA-"* ]]
            then
                allele=${all:0:7}":"${all:7:2}
            else
                allele=$all
            fi
            if [ $allele_count -lt 1 ]
            then
                cur_alleles=${cur_alleles}_${allele}
                allele_count=$((allele_count+1))
            else
                c=$((c+1))
                ln -sT ../peptides/$f $netmhc_dir/jobs/$c-$f.${cur_alleles:1} 2>/dev/null
                cur_alleles=_${allele}
                allele_count=1
            fi
        done
        if [ ! -z "$cur_alleles" ]
        then
            c=$((c+1))
            ln -sT ../peptides/$f $netmhc_dir/jobs/$c-$f.${cur_alleles:1} 2>/dev/null
        fi
    done
    popd
    
    # check for already completed jobs
    cc=
    shopt -s nullglob
    exist=(output/allPeptides_$genome/$netmhc_dir/done/*.done)
    shopt -u nullglob
    if [ -z "$exist" ]
    then
        ic=$c
        max=$((1000-njobs))
        if [ $ic -gt $max ]
        then
            ic=$max
        fi
        if [ $ic -gt 0 ]
        then
            cc=1-$ic
            njobs=$((njobs+ic))
        fi
    else
        for ic in `seq 1 $c` 
        do
            if [ $njobs -lt 1000 ]
            then
                if [ ! -f output/allPeptides_$genome/$netmhc_dir/done/$ic-*.done ]
                then
                    cc=$cc,$ic
                    njobs=$((njobs+1))
                fi
            fi
        done
        cc=`echo $cc | cut -c 2-`
    fi

    mkdir -p log/netmhcpan/tmp
    mkdir -p output/allPeptides_$genome/$netmhc_dir/predictions
    mkdir -p output/allPeptides_$genome/$netmhc_dir/done

    #echo $cc
    #echo $njobs
    
    if [ -z $cc ]
    then
        return
    fi

    echo '
    export TMPDIR=/tmp
    TMPDIR=`mktemp -d`
    mkdir -p $TMPDIR
    module load netmhcpan/4.0a
    file=`ls output/allPeptides_'$genome'/'$netmhc_dir'/jobs/${PBS_ARRAYID}-*`
    f=`basename $file`
    letter=`basename $file | cut -d "." -f 1 | cut -d "-" -f 2`
    length=`basename $file | cut -d "." -f 2`
    alleles=`basename $file | cut -d "." -f 3`
    for allele in ${alleles//_/ }
    do
        out=output/allPeptides_'$genome'/'$netmhc_dir'/predictions/$letter.$length.$allele.tsv
        tmpout=log/netmhcpan/tmp/$letter.$length.$allele.log
        netMHCpan -p -BA -a $allele -l $length -f $file -xls -xlsfile $out > $tmpout
    done
    touch output/allPeptides_'$genome'/'$netmhc_dir'/done/$f.done
    ' | qsub -l nodes=1:ppn=1,walltime=24:00:00,pmem=2gb -d $PWD -N $fullname -o log/netmhcpan -e log/netmhcpan -t $cc
}

runPBS GRCm38.98
runPBS GRCh38.98
