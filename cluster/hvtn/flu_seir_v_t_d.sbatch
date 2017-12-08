#!/bin/bash
#
# flu_vacc_v_t.sbatch
#
#SBATCH --job-name flu_seir_v_t_d  # Job name
#SBATCH --partition largemem     # Slurm partition to use
#SBATCH --ntasks 6               # Number of tasks to run (one CPU core will be allocated per task)
#SBATCH --time 10-11:59          # Wall time limit in D-HH:MM
#SBATCH --mem-per-cpu=8000       # Memory limit for each tasks (in MB)
#SBATCH --mail-type=END          # Type of email notification- BEGIN,END,FAIL,ALL
#SBATCH --mail-user=fintzij@uw.edu # Email to which notifications will be sent

# set outfile
#set outfile = "~/"flu_vacc_v_t_#.$$.Rout

# execute commands
module load R
R CMD BATCH --quiet --no-save --no-restore flu_seir_v_t_d.R flu_seir_v_t_d.$$.Rout
