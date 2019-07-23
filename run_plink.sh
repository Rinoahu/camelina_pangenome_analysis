#!/bin/bash
# use plink to perform gwas analysis
# this is only a demo

# merge all vcf, maybe need
bcftools merge --force-samples $1.vcf $2.vcf $3.vcf > merge.vcf


# convert vcf to files required by plink
plink --file merge.vcf --allow-no-sex --dog --make-bed --noweb --out merge.binary

# perform gwas
plink --bfile merge.binary --make-pheno merge.pheno "1" --assoc --allow-no-sex --adjust --dog --noweb --out merge_out

