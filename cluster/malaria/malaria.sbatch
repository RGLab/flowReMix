#!/bin/bash
#SBATCH --job-name mal_long # Set a name for your job. This is useful if you have multiple jobs queued.
#SBATCH --partition medium    # Slurm partition to use
#SBATCH --ntasks 16          # Number of tasks to run. By default, one CPU core will be allocated per task
#SBATCH --time 0-24:00       # Wall time limit in D-HH:MM
#SBATCH --mem-per-cpu=8000    # Memory limit for each tasks (in MB)
#SBATCH -o flowout.out     # File to which STDOUT will be written
#SBATCH -e flowout.out
#SBATCH --mail-type=ALL       # Type of email notification- BEGIN,END,FAIL,ALL
#SBATCH --mail-user=amitmeir@uw.edu # Email to which notifications will be sent

R CMD BATCH malariaCluster.R
