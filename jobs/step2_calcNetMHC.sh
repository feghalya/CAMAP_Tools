#!/usr/bin/env bash

name=calc_netMHC
date=`date +%Yy%mm%dd_%Hh%Mm%Ss`

# number of maximum jobs allowed
maxjobs=1000

# number of already running jobs
njobs=0

# do not submit job ids lower than this
min=1

# torque or slurm
if command -v qsub >/dev/null 2>&1
then
    scheduler_array_var=PBS_ARRAYID
else
    scheduler_array_var=SLURM_ARRAY_TASK_ID
fi


function runPBS {
    genome=$1
    NP=$2  # suffix to use if binding score (BA) not required

    if [ -z "$NP" ]; then BA='-BA '; else BA=; fi

    fullname=${name}_${genome}
    echo $fullname

    netmhc_dir=NetMHCpan-4.0a

    if [ -f output/allPeptides_$genome/$netmhc_dir.tar.gz ]
    then
        echo "Analysis already completed for genome $genome"
        return
    fi

    pushd output/allPeptides_$genome

    if [ ! -d $netmhc_dir/peptides ]
    then
        mkdir -p $netmhc_dir/peptides
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
    elif [[ $genome == "GRCm"* ]]
    then
        alleles=`cat ../../data/alleles/Mouse.tsv | cut -f 1 | tail -n +2 | head -n -1`
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
    exist=(output/allPeptides_$genome/$netmhc_dir/done/*[0-9bd]${NP}.done)
    shopt -u nullglob
    if [ -z "$exist" ]
    then
        ic=$c
        max=$((maxjobs+min-1-njobs))
        if [ $ic -gt $max ]
        then
            ic=$max
        fi
        if [ $ic -gt 0 ]
        then
            cc=$min-$ic
            njobs=$((njobs+ic-min+1))
        fi
    else
        for ic in `seq 1 $c` 
        do
            if [ $njobs -lt $maxjobs ]
            then
                if [ $ic -ge $min ] && [ ! -f output/allPeptides_$genome/$netmhc_dir/done/$ic-*[0-9bd]${NP}.done ]
                then
                    if [ -z "$cc" ]
                    then
                        cc=$ic-$ic
                    elif [ $(expr $(echo $cc | rev | cut -d '-' -f 1 | rev) + 1) -eq $ic ]
                    then
                        cc=$(echo $cc | rev | cut -d '-' -f 2- | rev)-$ic
                    else
                        cc=$cc,$ic-$ic
                    fi
                    njobs=$((njobs+1))
                fi
            fi
        done
    fi

    mkdir -p log/netmhcpan/tmp
    mkdir -p output/allPeptides_$genome/$netmhc_dir/predictions
    mkdir -p output/allPeptides_$genome/$netmhc_dir/done

    echo "jobs: $cc"
    echo "njobs: $njobs"

    if [ -z $cc ]
    then
        return
    fi

    cmd='#!/usr/bin/env bash''
    echo $'$scheduler_array_var'
    export TMPDIR=/tmp
    TMPDIR=`mktemp -d`
    mkdir -p $TMPDIR
    module load netmhcpan/4.0a
    file=`ls output/allPeptides_'$genome'/'$netmhc_dir'/jobs/${'$scheduler_array_var'}-*`
    f=`basename $file`
    letter=`basename $file | cut -d "." -f 1 | cut -d "-" -f 2`
    length=`basename $file | cut -d "." -f 2`
    alleles=`basename $file | cut -d "." -f 3`
    for allele in ${alleles//_/ }
    do
        out=output/allPeptides_'$genome'/'$netmhc_dir'/predictions/$letter.$length.$allele'$NP'.tsv
        tmpout=log/netmhcpan/tmp/$letter.$length.$allele'$NP'.log
        netMHCpan -p '$BA'-a $allele -l $length -f $file -xls -xlsfile $out > $tmpout
    done
    touch output/allPeptides_'$genome'/'$netmhc_dir'/done/$f'$NP'.done
    '

    echo "$cmd"

    if [ $scheduler_array_var = "PBS_ARRAYID" ]
    then
        echo "$cmd" |\
             qsub -l nodes=1:ppn=1:bioinfo,walltime=24:00:00,pmem=2gb -d $PWD -N $fullname \
                 -o log/netmhcpan -e log/netmhcpan -t $cc
    elif [ $scheduler_array_var = "SLURM_ARRAY_TASK_ID" ]
    then
        echo "$cmd" |\
            sbatch --export ALL --account $RAP_ID --chdir $PWD --time 0-12:00:00 --nodes 1 --cpus-per-task 1 \
                --mem-per-cpu 2gb --array $cc --output log/netmhcpan/$fullname.$date.%A_%j_%a.log \
                --error log/netmhcpan/$fullname.$date.%A_%j_%a.err --job-name $fullname
    fi
}


runPBS GRCh37.75
runPBS GRCh37.75 .NP
runPBS GRCm38.78
runPBS GRCm38.78 .NP
runPBS GRCh38.98
runPBS GRCm38.98
