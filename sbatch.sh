#!/bin/bash
#SBATCH --get-user-env
#SBATCH --job-name="RCA"
#SBATCH --time=10:00:00
#SBATCH --nodes=1
#SBATCH --mem=50GB
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
module load matlab/R2018b
module load perl/5.24.0

perl randomTF_target.pl