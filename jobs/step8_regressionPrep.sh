#!/usr/bin/env bash

name=regression_setup
workers=40
date=`date +%Yy%mm%dd_%Hh%Mm%Ss`

function runPBS {
    extra_params=$1
    scratch_location=$2
    params="${extra_params} -f ${scratch_location} "
    fullname=${name}${extra_params// /}
    echo $fullname

    echo '#!/usr/bin/env bash'"
module load scipy-stack/2021a
source /home/feghalya/.virtualenvs/ml2/bin/activate
time srun -n 1 ./regressionSetup.py -w $(($workers-1)) ${params}
sleep 5
sacct --format=JobID%15,State,ExitCode,CPUTime,MaxRSS,Start,End --units M -j \$SLURM_JOBID
" | sbatch --account $RAP_ID --chdir $PWD --time 0-1:00:00 --ntasks 1 --cpus-per-task $workers \
           --mem-per-cpu 18gb --output log/$fullname.$date.log --error log/$fullname.$date.err --job-name $fullname
}


module purge
mkdir -p log

runPBS "-n -r 1 -a SGD -p e4000 -m 10 -x 12 -s logregcv,rf,knn" "/scratch/feghalya/camap/f98aa5a2"
runPBS "-n -r 1 -a SGD -p e4000 -m 10 -x 12 -s mlp --balance" "/scratch/feghalya/camap/f98aa5a2"
runPBS "-n -r 1 -a SGD -p e4000 -m 10 -x 12 -s logregcv,rf,knn --max-tpm" "/scratch/feghalya/camap/e1fbb414"
runPBS "-n -r 1 -a SGD -p e4000 -m 10 -x 12 -s mlp --max-tpm --balance" "/scratch/feghalya/camap/e1fbb414"
