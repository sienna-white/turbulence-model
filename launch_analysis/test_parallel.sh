#!/bin/bash
#SBATCH --nodes 1
#SBATCH -N 1
#SBATCH -C cpu
#SBATCH -q regular
#SBATCH --ntasks-per-node=64
#SBATCH -J pmax_v_ws
#SBATCH -t 00:35:00
#SBATCH 

##### SBATCH -account m1266

module load PrgEnv-gnu
module load python

cd .. 

python watercolumn-parallel.py
