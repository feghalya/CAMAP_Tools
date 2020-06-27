#!/usr/bin/env bash

name=merge_netMHC
date=`date +%Yy%mm%dd_%Hh%Mm%Ss`

ntasks=40

function runPBS {
    genome=$1
    fullname=${name}_${genome}
    echo $fullname

    echo '#!/usr/bin/env bash'"
time srun -n $ntasks python -m mpi4py.futures ./mergeNetMHC.py -w $(($ntasks-1)) -g $genome
sleep 5
sacct --format=JobID%15,State,ExitCode,CPUTime,MaxRSS,Start,End --units M -j \$SLURM_JOBID
" | sbatch --account $RAP_ID --workdir $PWD --time 0-12:00:00 --ntasks $ntasks --cpus-per-task 1 \
	   --mem-per-cpu 32gb --output log/$fullname.$date.log --error log/$fullname.$date.err --job-name $fullname

}

#runPBS GRCh38.98
runPBS GRCm38.98
