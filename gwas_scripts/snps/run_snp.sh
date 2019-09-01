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
#SBATCH --mem         128000         # 2000 MB of Memory allocated; set --mem with care
#SBATCH --time        24:0:0     # Maximum job run time
##SBATCH --mail-user   $email # user to send emails to
##SBATCH --mail-type   ALL                 # Email on: BEGIN, END, FAIL & REQUEUE
## Run 'man sbatch' for more information on the options above.


cd /mnt/lustrefs/scratch/xiao.hu1/camelina/JGI_full/bams/snps

#bcftools mpileup -m 3 -F 0.0002 -uf $ref $bams > out.bcf
# | bcftools view -bvcg - > out.bcf
#echo $bams
#bcftools mpileup --threads 12 -Ou -f $ref $bams > call.vcf

ref="../../tests/Cs_genome_v2.fa"
/usr/bin/time -v bcftools mpileup --threads 6 -Ou -f $ref -b ./list.txt | bcftools call --threads 6 -mv -Ob -o calls_all.bcf &> log_all.txt

