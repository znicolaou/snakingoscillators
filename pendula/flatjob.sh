#!/bin/bash
#SBATCH --account=amath
#SBATCH --partition=ckpt
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=20G
#SBATCH --time=01-00:00:00 # Max runtime in DD-HH:MM:SS format.
#SBATCH --export=all
#SBATCH --output=flat.out 
#SBATCH --error=flat.err 

mkdir -p data/pendula/flat/
./makeflat.py --num 16 --file data/pendula/flat/flat.dat
cp c.flat data/pendula/flat/
cp pendula.c data/pendula/flat/
cp flat.auto data/pendula/flat/
cd data/pendula/flat
auto flat.auto
cd ..
for i in {0..5}; do
	mv flat/${i} ./$((i+19))
done
