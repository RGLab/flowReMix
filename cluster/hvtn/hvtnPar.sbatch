#!/bin/bash
#SBATCH --job-name HVTN # Set a name for your job. This is useful if you have multiple jobs queued.
#SBATCH --partition short    # Slurm partition to use
#SBATCH --ntasks 9         # Number of tasks to run. By default, one CPU core will be allocated per task
#SBATCH --time 0-12:00        # Wall time limit in D-HH:MM
#SBATCH --mem-per-cpu=7000    # Memory limit for each tasks (in MB)
#SBATCH -o out/flowout_%j.out
#SBATCH -e err/flowout_%j.out
#SBATCH --mail-type=ALL       # Type of email notification- BEGIN,END,FAIL,ALL
#SBATCH --mail-user=amitmeir@uw.edu # Email to which notifications will be sent
#SBATCH --array=1-400

module load R

R CMD BATCH --no-save --no-restore "--args setting='${SLURM_ARRAY_TASK_ID}'" HVTNclusterSA7.R

