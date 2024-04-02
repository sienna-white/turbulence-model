#!/bin/sh
#BATCH --job-name=run_analysis
#SBATCH --partition=savio3 
##SBATCH --qos=aiolos_savio3_normal 
#SBATCH --account=co_aiolos
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
#### // SBATCH --cpus-per-task=20
#SBATCH --time=00:20:59


# module load PrgEnv-gnu
module load python




cd .. 




python watercolumn-parallel.py
