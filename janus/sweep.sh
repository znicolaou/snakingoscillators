#!/bin/bash
#SBATCH --account=amath
#SBATCH --partition=ckpt
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=20G
#SBATCH --time=01-00:00:00 # Max runtime in DD-HH:MM:SS format.
#SBATCH --export=all
#SBATCH --array=0-100
#SBATCH --output=out/%a.out 
#SBATCH --error=out/%a.err 

seeds=100
seed0=$((seeds*SLURM_ARRAY_TASK_ID))
export OMP_NUM_THREADS=1;

for i in `seq 1 $seeds`; do 
	seed=$((seed0+i));
	jobs=`jobs | wc -l`;
	while [ $jobs -ge 16 ]; do
		sleep 1;
		jobs=`jobs | wc -l`;
	done
	if [ ! -f data/randomjanus/${seed}out.dat ]; then
		echo $seed
		./janus.py --filebase data/randomjanus/$seed --seed $seed --sigma 0.3 --beta 0.3 --dt 0.01 --sym 0 --num 16 --time 10000 --rtime 5000 --output 0 &
	fi
done
wait

