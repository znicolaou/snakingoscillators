#!/bin/bash
#SBATCH --account=amath
#SBATCH --partition=ckpt
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=5G
#SBATCH --time=02-00:00:00 # Max runtime in DD-HH:MM:SS format.
#SBATCH --export=all
#SBATCH --output=%a/out.dat # where STDOUT goes
#SBATCH --error=%a/err.dat # where STDERR goes

# Your programs to run.
cd $SLURM_ARRAY_TASK_ID
cp ../cycles.auto ./
cp ../c.forward ./
cp ../janus_rel.c ./
auto cycles.auto
