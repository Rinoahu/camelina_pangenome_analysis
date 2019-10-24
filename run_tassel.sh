#!/bin/bash
##
## hello.slurm.sh: a simple slurm batch job
##
## Lines starting with #SBATCH are read by Slurm. Lines starting with ## are comments.
## All other lines are read by the shell.
##
#SBATCH --job-name    tassel        # job name
#SBATCH --output      tassel-%j.out # standard output file (%j = jobid)
#SBATCH --error       tassel-%j.err # standard error file
##SBATCH --partition   express      # queue partition to run the job in
#SBATCH --nodes       1            # number of nodes to allocate
#SBATCH --ntasks-per-node 32        # number of cores to allocate; set with care 
#SBATCH --mem         128000         # 2000 MB of Memory allocated; set --mem with care
#SBATCH --time        18:0:0     # Maximum job run time
##SBATCH --mail-user   $email # user to send emails to
##SBATCH --mail-type   ALL                 # Email on: BEGIN, END, FAIL & REQUEUE
## Run 'man sbatch' for more information on the options above.


cd /mnt/lustrefs/scratch/xiao.hu1/camelina/JGI_full/tassel_output


#vcftools --vcf merge_all.vcf --max-missing 0.2 --mac 3 --minQ 30 --remove-indels --out merge_all.snps --recode --recode-INFO-all

beagle nthreads=32 gt=merge_all.snps.recode.vcf out=merge_all.snps.impute_bg chrom=Chr17

#plink --vcf merge_all.snps.impute_bg.vcf.gz --maf .05 --geno .5 --recode vcf-iid --out merge_all.snps.impute_bg_flt --allow-extra-chr
#plink --vcf merge_all.snps.impute_bg_flt --indep-pairwise 50 10 .2 --out test_filter_LD --allow-extra-chr


#run_pipeline.pl -Xms64G -Xmx125g -importGuess merge_all.vcf -LDKNNiImputationPlugin -highLDSSites 50 -knnTaxa 10 -maxLDDistance 100000000 -endPlugin -export merge_all0.imputed.vcf -exportType VCF -sortPositions

#run_pipeline.pl -Xms64G -Xmx125g -importGuess merge_0.vcf -LDKNNiImputationPlugin -highLDSSites 50 -knnTaxa 10 -maxLDDistance 100000000 -endPlugin -export merge_0.imputed.vcf -exportType VCF -sortPositions


#beagle nthreads=32 gt=merge_all.vcf out=merge_all_impute_bg chrom=Chr17

#plink --vcf merge_0.imputed.vcf --maf .05 --geno .5 --recode vcf-iid --out test_filter0 --allow-extra-chr
#plink --vcf test_filter0.vcf --indep-pairwise 50 10 .2 --out test_filter_LD --allow-extra-chr


#run_pipeline.pl -Xms10g -Xmx100g  -vcf Camsat_diversity_filtered_SNPs.vcf.gz -sortPositions -export out.hmp.txt -exportType HapmapDiploid
