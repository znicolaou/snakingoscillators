#!/bin/bash
#SBATCH --account=amath
#SBATCH --partition=ckpt
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=5G
#SBATCH --time=02-00:00:00 # Max runtime in DD-HH:MM:SS format.
#SBATCH --export=all
#SBATCH --output=data/pendula/%a/out.dat # where STDOUT goes
#SBATCH --error=data/pendula/%a/err.dat # where STDERR goes
#SBATCH --array=0-18

# Your programs to run.
cd data/pendula/$SLURM_ARRAY_TASK_ID
cp ../../../cycles.auto ./
cp ../../../c.forward ./
cp ../../../pendula.c ./
auto cycles.auto
grep 'Time' out.dat | awk '{s=s+$3}END{print("Runtime:", s)}'
