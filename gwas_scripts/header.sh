#!/bin/bash
##
## hello.slurm.sh: a simple slurm batch job
##
## Lines starting with #SBATCH are read by Slurm. Lines starting with ## are comments.
## All other lines are read by the shell.
##
#SBATCH --job-name    filter        # job name
#SBATCH --output      filter-%j.out # standard output file (%j = jobid)
#SBATCH --error       filter-%j.err # standard error file
##SBATCH --partition   express      # queue partition to run the job in
#SBATCH --nodes       1            # number of nodes to allocate
#SBATCH --ntasks-per-node 24        # number of cores to allocate; set with care 
#SBATCH --mem         96000         # 2000 MB of Memory allocated; set --mem with care
#SBATCH --time        48:0:0     # Maximum job run time
##SBATCH --mail-user   $email # user to send emails to
##SBATCH --mail-type   ALL                 # Email on: BEGIN, END, FAIL & REQUEUE
## Run 'man sbatch' for more information on the options above.

