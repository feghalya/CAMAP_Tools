#!/usr/bin/env bash

name=annotate_proteome
workers=24
date=`date +%Yy%mm%dd_%Hh%Mm%Ss`

function runPBS {
    genome=$1
    fullname=${name}_${genome}
    echo $fullname
    filters=$2
    if [ ! -z $filters ]; then filters="-f $filters "; fi

    echo '#!/usr/bin/env bash'"
source /home/feghalya/.virtualenvs/ml/bin/activate
time srun -n $workers python -m mpi4py.futures ./annotateProteome.py -g $genome -w $(($workers-1)) ${filters}--mpi
sleep 5
sacct --format=JobID%15,State,ExitCode,CPUTime,MaxRSS,Start,End --units M -j \$SLURM_JOBID
" | sbatch --account $RAP_ID --chdir $PWD --time 0-6:00:00 --ntasks $workers --cpus-per-task 8 \
           --mem-per-cpu 4gb --output log/$fullname.$date.log --error log/$fullname.$date.err --job-name $fullname
}


module purge
runPBS GRCh37.75 Adam,Adam_Shuffle
runPBS GRCm38.78 Adam,Adam_Shuffle
#runPBS GRCh37.75 SGD,SGD_Shuffle
#runPBS GRCm38.78 SGD,SGD_Shuffle
#runPBS GRCh38.98 Adam,Adam_Shuffle
#runPBS GRCm38.98 Adam,Adam_Shuffle
