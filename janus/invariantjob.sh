#!/bin/bash
#SBATCH --account=amath
#SBATCH --partition=ckpt
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=5G
#SBATCH --time=02-00:00:00 # Max runtime in DD-HH:MM:SS format.
#SBATCH --export=all
#SBATCH --output=out.dat # where STDOUT goes
#SBATCH --error=err.dat # where STDERR goes

i=8
br='rv5'
filebase='data/janus'
mkdir -p $filebase/$i/invariant
cp {janus_rel.c,c.forward,invariant.auto} $filebase/$i/invariant
cp $filebase/$i/s.start_$br $filebase/$i/invariant/s.start
cp $filebase/$i/b.start_$br $filebase/$i/invariant/b.start
cp $filebase/$i/d.start_$br $filebase/$i/invariant/d.start
cd $filebase/$i/invariant 
auto invariant.auto
