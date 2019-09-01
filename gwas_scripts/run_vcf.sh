#!/bin/bash
##
## hello.slurm.sh: a simple slurm batch job
##
## Lines starting with #SBATCH are read by Slurm. Lines starting with ## are comments.
## All other lines are read by the shell.
##
#SBATCH --job-name    merge        # job name
#SBATCH --output      merge-%j.out # standard output file (%j = jobid)
#SBATCH --error       merge-%j.err # standard error file
##SBATCH --partition   express      # queue partition to run the job in
#SBATCH --nodes       1            # number of nodes to allocate
#SBATCH --ntasks-per-node 12        # number of cores to allocate; set with care 
#SBATCH --mem         16000         # 2000 MB of Memory allocated; set --mem with care
#SBATCH --time        24:0:0     # Maximum job run time
##SBATCH --mail-user   $email # user to send emails to
##SBATCH --mail-type   ALL                 # Email on: BEGIN, END, FAIL & REQUEUE
## Run 'man sbatch' for more information on the options above.

#cd /mnt/lustrefs/scratch/xiao.hu1/camelina/JGI_full/vcfs
cd /mnt/lustrefs/scratch/xiao.hu1/camelina/JGI_full/tests
#echo $PWD > hello.txt
#echo $idba_hybrid >> hello.txt

x=''
y=''
for i in `ls *.vcf.gz`
do
	#bcftools index $i
	x="$x $i"
	y="$y -V $i"
done
#echo $y

#/usr/bin/time -v bcftools merge -m all --threads 12 --force-samples -f PASS $x -o merge.vcf -Oz &> log2.txt

/usr/bin/time -v $gatk3 -T CombineVariants -R ./Cs_genome_v2.fa $y -o output.vcf -genotypeMergeOptions UNIQUIFY &>log.txt

#python chn.py merge2.vcf > merge2_cn.vcf
