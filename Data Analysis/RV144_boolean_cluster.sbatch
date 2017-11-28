#!/bin/bash
#SBATCH --job-name RV144 # Set a name for your job. This is useful if you have multiple jobs queued.
#SBATCH --partition=largenode    # Slurm partition to use
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --tasks-per-node=1 
#SBATCH --mem=33000
#SBATCH --time 0-24:00        # Wall time limit in D-HH:MM
#SBATCH -o logs/rv144_%A_%a.out     # File to which STDOUT will be written
#SBATCH -e errs/rv144_err_%A_%a.out
#SBATCH --mail-type=ALL       # Type of email notification- BEGIN,END,FAIL,ALL
#SBATCH --mail-user=gfinak@fredhutch.org # Email to which notifications will be sent
#SBATCH --array=1-60
#SBATCH --cpus-per-task=16

Rscript --no-save --no-restore RV144_boolean_cluster.R 
