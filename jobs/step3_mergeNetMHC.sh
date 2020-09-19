#!/usr/bin/env bash

name=merge_netMHC
date=`date +%Yy%mm%dd_%Hh%Mm%Ss`

ntasks=8

function runPBS {
    genome=$1
    fullname=${name}_${genome}
    echo $fullname

    echo '#!/usr/bin/env bash'"
source /home/feghalya/.virtualenvs/ml/bin/activate
time srun -n $ntasks python -m mpi4py.futures ./mergeNetMHC.py -w $(($ntasks-1)) -g $genome --mpi
sleep 5
sacct --format=JobID%15,State,ExitCode,CPUTime,MaxRSS,Start,End --units M -j \$SLURM_JOBID
" | sbatch --account $RAP_ID --chdir $PWD --time 0-12:00:00 --ntasks $ntasks --cpus-per-task 1 \
	   --mem-per-cpu 32gb --output log/$fullname.$date.log --error log/$fullname.$date.err --job-name $fullname

}

module purge

runPBS GRCh37.75
#runPBS GRCm38.78
#runPBS GRCh38.98
#runPBS GRCm38.98
