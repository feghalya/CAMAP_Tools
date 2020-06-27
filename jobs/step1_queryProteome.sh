#!/usr/bin/env bash

name=query_proteome
workers=16
date=`date +%Yy%mm%dd_%Hh%Mm%Ss`

function runPBS {
    genome=$1
    fullname=${name}_${genome}
    echo $fullname

    echo '#!/usr/bin/env bash'"
cp -rT /home/feghalya/.pyGeno.bak /dev/shm/pyGeno
source /home/feghalya/.virtualenvs/ml/bin/activate
time srun -n 1 ./queryAllProteome.py -c 162 -w $(($workers-1)) -g $genome
sleep 5
rm -fr /dev/shm/pyGeno
sacct --format=JobID%15,State,ExitCode,CPUTime,MaxRSS,Start,End --units M -j \$SLURM_JOBID
" | sbatch --export ALL --account $RAP_ID --workdir $PWD --time 0-6:00:00 --nodes 1 --cpus-per-task $workers \
	   --mem-per-cpu 16gb --output log/$fullname.$date.log --error log/$fullname.$date.err --job-name $fullname
}

mkdir -p log

runPBS GRCh38.98
runPBS GRCm38.98
