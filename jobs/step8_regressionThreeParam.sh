#!/usr/bin/env bash

name=regression_analysis
workers=25
date=`date +%Yy%mm%dd_%Hh%Mm%Ss`

function runPBS {
    extra_params=$1
    fullname=${name}${extra_params// /}
    echo $fullname

    echo '#!/usr/bin/env bash'"
source /home/feghalya/.virtualenvs/ml/bin/activate
time srun -n $workers python -m mpi4py.futures ./regressionAnalysis.py -w $(($workers-1)) ${extra_params}--mpi
sleep 5
sacct --format=JobID%15,State,ExitCode,CPUTime,MaxRSS,Start,End --units M -j \$SLURM_JOBID
" | sbatch --account $RAP_ID --chdir $PWD --time 0-6:00:00 --ntasks $workers --cpus-per-task 16 \
           --mem-per-cpu 4gb --output log/$fullname.$date.log --error log/$fullname.$date.err --job-name $fullname
}


module purge
mkdir -p log

#runPBS "-r 1250 -a Adam -p e500pep9t5 "
runPBS "-r 1 -a Adam -p e500pep9t5 "
#runPBS "-n -r 1 -a Adam -p e500pep9t5 "
