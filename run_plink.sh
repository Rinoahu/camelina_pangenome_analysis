#!/bin/bash
# use plink to perform gwas analysis
# this is only a demo

# merge all vcf, maybe need
bcftools merge --force-samples $1.vcf $2.vcf $3.vcf > merge.vcf

bcftool query -l  merge.vcf > order_accession.txt


# convert vcf to files required by plink
plink --file merge.vcf --allow-no-sex --dog --make-bed --noweb --out merge.binary --allow-extra-chr

# the first 2 col of *.fam is the name for pheno

# perform gwas
plink --bfile merge.binary --make-pheno merge.pheno "1" --assoc --allow-no-sex --adjust --dog --noweb --out merge_out

