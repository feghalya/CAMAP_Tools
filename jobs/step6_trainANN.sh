#!/usr/bin/env bash

date=`date +%Yy%mm%dd_%Hh%Mm%Ss`

function runPBS {
    name=$1
    exp=$2
    dsname=$3
    run=$4

    echo $name
    echo $exp
    echo $dsname
    echo $run

    outdir=log/camap

    echo '#!/usr/bin/env bash'"
source /home/feghalya/.virtualenvs/ml/bin/activate
for s in {$exp..$((exp+2))}
do
    infile=output/trainDS/${dsname}_seed\$s
    echo '$run' -ds \$infile
    $run -ds \$infile 2>&1 > $outdir/$name-\$s-$date.log | tr '\r' '\n' > $outdir/$name-\$s-$date.err &
    sleep 5
done
wait
" | sbatch --export ALL --account $BLUES_ID --chdir $PWD --time 0-12:00:00 --nodes 1 --cpus-per-task 16 --mem-per-cpu 1gb --gres gpu:1 --output $outdir/$name$exp-$date.log --error $outdir/$name$exp-$date.err --job-name $name-$exp
}


module purge
mkdir -p log/camap

context=162
epochs=500

datasets=
#datasets=$datasets" BLCL_GRCh37.75_padding162_maxBS500_maxContexts10_ratio5_peplen9@pep9t5"
datasets=$datasets" BLCL_GRCh37.75_padding162_maxBS500_maxContexts10_ratio5_peplen9_sameTPM@pep9t5y"
datasets=$datasets" BLCL_GRCh37.75_padding162_maxBS500_maxContexts3_ratio5_peplen9@pep9t5m3"
datasets=$datasets" BLCL_GRCh37.75_padding162_maxBS500_maxContexts0_ratio5_peplen9@pep9t5mINF"
#datasets=$datasets" BLCL_GRCh38.98_padding162_maxBS500_maxContexts10_ratio5_peplen9@pep9t5"
#datasets=$datasets" BLCL_GRCh38.98_padding162_maxBS500_maxContexts10_ratio5_peplen9_sameTPM@pep9t5y"
#datasets=$datasets" BLCL_GRCh38.98_padding162_maxBS500_maxContexts10_ratio25_peplen9@pep9t25"
#datasets=$datasets" BLCL_GRCh38.98_padding162_maxBS500_maxContexts10_ratio25_peplen9_sameTPM@pep9t25y"
#datasets=$datasets" BLCL_GRCh38.98_padding162_maxBS500_maxContexts10_ratio5@pepallt5"

for dsname in $datasets
do
    dsshortname=${dsname#*@}
    dsname=${dsname%@*}

    for exp in 1 4 7 10
    do
        name="p${context}e${epochs}${dsshortname}"
        run="camap $context -e $epochs -n $name -w 5 -o adam --device cuda -subf pytorch-adam"
        runPBS "$name-adam" "$exp" "$dsname" "$run"
    done

    for exp in 1 4 7 10
    do
        name="p${context}e${epochs}${dsshortname}"
        run="camap ${context} -cse -e $epochs -n $name -w 5 -o adam --device cuda -subf pytorch-adam-shuffle"
        runPBS "$name-adam-shuffle" "$exp" "$dsname" "$run"
    done
done
