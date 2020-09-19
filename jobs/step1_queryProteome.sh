#!/usr/bin/env bash

name=query_proteome
workers=16
date=`date +%Yy%mm%dd_%Hh%Mm%Ss`

function runPBS {
    genome=$1
    fullname=${name}_${genome}
    echo $fullname

    echo '#!/usr/bin/env bash'"
set -e
stagein() { cp -rT /home/feghalya/.pyGeno.bak /dev/shm/pyGeno; ln -fsT /dev/shm/pyGeno /home/feghalya/.pyGeno; }
stageout() { rm -fr /dev/shm/pyGeno; ln -fsT /home/feghalya/.pyGeno.bak /home/feghalya/.pyGeno; echo 'Done!' >&2; }
trap stageout EXIT
stagein
source /home/feghalya/.virtualenvs/ml/bin/activate
time srun -n 1 ./queryAllProteome.py -c 162 -w $(($workers-1)) -g $genome
sleep 5
sacct --format=JobID%15,State,ExitCode,CPUTime,MaxRSS,Start,End --units M -j \$SLURM_JOBID
" | sbatch --export ALL --account $RAP_ID --chdir $PWD --time 0-6:00:00 --nodes 1 --cpus-per-task $workers \
	   --mem-per-cpu 16gb --output log/$fullname.$date.log --error log/$fullname.$date.err --job-name $fullname
}

module purge
mkdir -p log

runPBS GRCh37.75
runPBS GRCm38.78
runPBS GRCh38.98
runPBS GRCm38.98
