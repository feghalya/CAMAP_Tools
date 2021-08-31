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
module load mpi4py/3.0.3
source /home/feghalya/.virtualenvs/ml2/bin/activate
time srun -n $workers python -m mpi4py.futures ./regressionSetup.py -w $(($workers-1)) ${params}--mpi
sleep 5
sacct --format=JobID%15,State,ExitCode,CPUTime,MaxRSS,Start,End --units M -j \$SLURM_JOBID
" | sbatch --account $RAP_ID --chdir $PWD --time 0-6:00:00 --ntasks $workers --cpus-per-task 1 \
           --mem-per-cpu 24gb --output log/$fullname.$date.log --error log/$fullname.$date.err --job-name $fullname
}


module purge
mkdir -p log

runPBS "-n -r 1 -a SGD -p e4000 -m 10 -x 12 -s mlp,logregcv,rf,knn" "/scratch/feghalya/camap/39cf757a"
