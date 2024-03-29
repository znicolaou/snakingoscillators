#!/bin/bash
#SBATCH --account=amath
#SBATCH --partition=ckpt
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=5G
#SBATCH --time=02-00:00:00 # Max runtime in DD-HH:MM:SS format.
#SBATCH --export=all
#SBATCH --output=data/janus/%a/out.dat # where STDOUT goes
#SBATCH --error=data/janus/%a/err.dat # where STDERR goes
#SBATCH --array=1-64

# Your programs to run.
cd data/janus/$SLURM_ARRAY_TASK_ID
cp ../../../cycles.auto ./
cp ../../../c.forward ./
cp ../../../janus_rel.c ./
cp ../../../maketotal.auto ./
auto cycles.auto
#auto ../maketotal.auto
