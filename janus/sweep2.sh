#!/bin/bash
#SBATCH --account=amath
#SBATCH --partition=ckpt
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=20G
#SBATCH --time=01-00:00:00 # Max runtime in DD-HH:MM:SS format.
#SBATCH --export=all
#SBATCH --array=1-100
#SBATCH --output=/mmfs1/home/zgn/snakingoscillators/out/%a.out 
#SBATCH --error=/mmfs1/home/zgn/snakingoscillators/out/%a.err 

export OMP_NUM_THREADS=1;
if [ $# -ne 1 ]; then echo "need seeds"; exit; fi
seeds="$1"
mkdir -p data/randomjanus/cycles
for seed in $seeds; do 
	jobs=`jobs | wc -l`;
	while [ $jobs -ge 16 ]; do
		sleep 1;
		jobs=`jobs | wc -l`;
	done
	if [ ! -f data/randomjanus/cycles/${seed}cycle.dat ]; then
		cp data/randomjanus/${seed}fs.npy data/randomjanus/cycles/${seed}ic.npy
		echo $seed
		./janus.py --filebase data/randomjanus/cycles/$seed --seed $seed --sigma 0.3 --beta 0.3 --dt 0.01 --sym 0 --num 16 --time 2000 --rtime 0 --output 3 &
	fi
done
wait

