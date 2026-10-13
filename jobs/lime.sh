#!/usr/bin/env bash

name=lime
workers=12
date=`date +%Yy%mm%dd_%Hh%Mm%Ss`

function runPBS {
    echo '#!/usr/bin/env bash'"
source /home/feghalya/.virtualenvs/ml/bin/activate
time srun -n $workers python -m mpi4py.futures ./utils/lime_prefs.py
sleep 5
sacct --format=JobID%15,State,ExitCode,CPUTime,MaxRSS,Start,End --units M -j \$SLURM_JOBID
" | sbatch --account $RAP_ID --chdir $PWD --time 0-24:00:00 --ntasks $workers --cpus-per-task 16 \
           --mem-per-cpu 4gb --output $name.$date.log --error $name.$date.err --job-name $name
}

module purge
runPBS
