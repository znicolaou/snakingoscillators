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
#SBATCH --output=/mmfs1/home/zgn/snakingoscillators/out/%a.out 
#SBATCH --error=/mmfs1/home/zgn/snakingoscillators/out/%a.err 

export OMP_NUM_THREADS=1;

for seed in $@; do 
	jobs=`jobs | wc -l`;
	while [ $jobs -ge 16 ]; do
		sleep 1;
		jobs=`jobs | wc -l`;
	done
	if [ ! -f data/randompendula/${seed}dat.npy ]; then
		echo $seed
		./pendula.py --frequency 3.5 --amplitude 0.045 --delta 0.25 --init 0.5 --cycles 5000 --outcycle 4900 --dt 0.01 --num 32 --seed $seed --filebase data/randompendula/$seed --verbose 1 &
	fi
done
wait

