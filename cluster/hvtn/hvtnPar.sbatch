#!/bin/bash
#SBATCH --job-name HVTN
#SBATCH --partition short
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --time 0-12:00        # Wall time limit in D-HH:MM
#SBATCH --mem-per-cpu=4000    # Memory limit for each tasks (in MB)
#SBATCH -o out/flowout_%j.out
#SBATCH -e err/flowout_%j.out
#SBATCH --mail-type=ALL       # Type of email notification- BEGIN,END,FAIL,ALL
#SBATCH --mail-user=amitmeir@uw.edu # Email to which notifications will be sent
#SBATCH --array=1-80

module load R

R CMD BATCH --no-save --no-restore "--args setting='${SLURM_ARRAY_TASK_ID}'" HVTNclusterSA8.R

