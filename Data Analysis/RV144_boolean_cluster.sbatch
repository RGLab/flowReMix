#!/bin/bash
#SBATCH --job-name RV144 # Set a name for your job. This is useful if you have multiple jobs queued.
#SBATCH --partition medium    # Slurm partition to use
#SBATCH --ntasks 10         # Number of tasks to run. By default, one CPU core will be allocated per task
#SBATCH --time 0-24:00        # Wall time limit in D-HH:MM
#SBATCH --mem-per-cpu=4500    # Memory limit for each tasks (in MB)
#SBATCH -o rv144.out     # File to which STDOUT will be written
#SBATCH -e rv144.out
#SBATCH --mail-type=ALL       # Type of email notification- BEGIN,END,FAIL,ALL
#SBATCH --mail-user=gfinak@fredhutch.org # Email to which notifications will be sent
#SBATCH --array=1-60
#SBATCH --cpus-per-task=20

module load R

R CMD BATCH --no-save --no-restore RV144_boolean_cluster.R
