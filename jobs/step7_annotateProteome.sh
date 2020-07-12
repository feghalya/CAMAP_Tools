#!/usr/bin/env bash

name=annotate_proteome
workers=24
date=`date +%Yy%mm%dd_%Hh%Mm%Ss`

function runPBS {
    genome=$1
    name=${name}_${genome}
    echo $name

    echo '#!/usr/bin/env bash'"
source /home/feghalya/.virtualenvs/ml/bin/activate
time srun -n $workers python -m mpi4py.futures ./annotateProteome.py -g $genome -w $(($workers-1)) --mpi
sleep 5
sacct --format=JobID%15,State,ExitCode,CPUTime,MaxRSS,Start,End --units M -j \$SLURM_JOBID
" | sbatch --account $RAP_ID --workdir $PWD --time 0-12:00:00 --ntasks $workers \
           --mem-per-cpu 32gb --output log/$name.$date.log --error log/$name.$date.err --job-name $name
}


module purge
runPBS GRCh37.75
runPBS GRCm38.78
