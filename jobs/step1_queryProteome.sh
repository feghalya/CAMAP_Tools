#!/usr/bin/env bash

name=query_proteome
workers=16
date=`date +%Yy%mm%dd_%Hh%Mm%Ss`

function runPBS {
    # assumes ~/.pyGeno is a link to /dev/shm/pyGeno
    genome=$1
    name=${name}_${genome}
    echo $name

    echo '#!/usr/bin/env bash'"
source /home/feghalya/.virtualenvs/ml/bin/activate
time srun -n 1 ./queryAllProteome.py -c 162 -w $(($workers-1)) -g $genome
sleep 5
sacct --format=JobID%15,State,ExitCode,CPUTime,MaxRSS,Start,End --units M -j \$SLURM_JOBID
" | sbatch --export ALL --account $RAP_ID --workdir $PWD --time 0-12:00:00 --nodes 1 --cpus-per-task $workers \
	   --mem-per-cpu 16gb --output log/$name.$date.log --error log/$name.$date.err --job-name $name
}

runPBS GRCh37.75
runPBS GRCm38.78
