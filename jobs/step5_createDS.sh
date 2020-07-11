#!/usr/bin/env bash

name=create_ds
workers=24
date=`date +%Yy%mm%dd_%Hh%Mm%Ss`

function runPBS {
    cmd=$1
    n=$2
    fullname=${name}_${n}
    echo $fullname

    echo '#!/usr/bin/env bash'"
source /home/feghalya/.virtualenvs/ml/bin/activate
echo '$cmd'
time srun -n 1 $cmd
sleep 5
sacct --format=JobID%15,State,ExitCode,CPUTime,MaxRSS,Start,End --units M -j \$SLURM_JOBID
" | sbatch --export ALL --account $RAP_ID --workdir $PWD --time 0-2:00:00 --nodes 1 --cpus-per-task $workers \
	   --mem-per-cpu 12gb --output log/$fullname.$date.log --error log/$fullname.$date.err --job-name $fullname
}

module purge
mkdir -p log

genome=GRCh38.98
dataset=BLCL

job_cmd="./constructTrainDS.py -w $(($workers-1)) -g $genome -d $dataset -s 12 -c 162 -r 500 -m 10 -t 5 -p 9"
job_name=${genome}_${dataset}_p9_t5
runPBS "$job_cmd" "$job_name"

job_cmd="./constructTrainDS.py -w $(($workers-1)) -g $genome -d $dataset -s 12 -c 162 -r 500 -m 10 -t 5 -p 9 -y"
job_name=${genome}_${dataset}_p9_t5_y
runPBS "$job_cmd" "$job_name"

job_cmd="./constructTrainDS.py -w $(($workers-1)) -g $genome -d $dataset -s 12 -c 162 -r 500 -m 10 -t 25 -p 9"
job_name=${genome}_${dataset}_p9_t25
runPBS "$job_cmd" "$job_name"

job_cmd="./constructTrainDS.py -w $(($workers-1)) -g $genome -d $dataset -s 12 -c 162 -r 500 -m 10 -t 25 -p 9 -y"
job_name=${genome}_${dataset}_p9_t25_y
runPBS "$job_cmd" "$job_name"

job_cmd="./constructTrainDS.py -w $(($workers-1)) -g $genome -d $dataset -s 12 -c 162 -r 500 -m 10 -t 5"
job_name=${genome}_${dataset}_pall_t5
runPBS "$job_cmd" "$job_name"
