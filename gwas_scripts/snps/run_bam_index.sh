#!/bin/bash
##
## hello.slurm.sh: a simple slurm batch job
##
## Lines starting with #SBATCH are read by Slurm. Lines starting with ## are comments.
## All other lines are read by the shell.
##
#SBATCH --job-name    snp        # job name
#SBATCH --output      snp-%j.out # standard output file (%j = jobid)
#SBATCH --error       snp-%j.err # standard error file
##SBATCH --partition   express      # queue partition to run the job in
#SBATCH --nodes       1            # number of nodes to allocate
#SBATCH --ntasks-per-node 12        # number of cores to allocate; set with care 
#SBATCH --mem         64000         # 2000 MB of Memory allocated; set --mem with care
#SBATCH --time        16:0:0     # Maximum job run time
##SBATCH --mail-user   $email # user to send emails to
##SBATCH --mail-type   ALL                 # Email on: BEGIN, END, FAIL & REQUEUE
## Run 'man sbatch' for more information on the options above.

module load parallel/20150622

cd /mnt/lustrefs/scratch/xiao.hu1/camelina/JGI_full/bams/snps

rm run_parallel_index.sh

for i in `find ../ -name *.bam | grep -v map | grep Resequencing`
do
	echo samtools index $i >> run_parallel_index.sh
done

parallel -j 16 < run_parallel_index.sh

