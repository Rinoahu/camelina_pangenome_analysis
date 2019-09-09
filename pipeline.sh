#!/bin/bash

########################################################################################
# survey of the first sequenced genome of camelina
########################################################################################

# download genome data
wget -c http://www.camelinadb.ca/downloads/Cs_genome_v2.fa.gz
tar xvf Cs_genome_v2.fa.gz


# kmer statistics
m=31
jellyfish count -m $m -s 200M -t 4 -o $1\.bin $1
jellyfish dump $1\.bin -o $1\.bin.txt



#########################################################################################
# resequencing
#########################################################################################

# pipeline for de novo assembly
# 1. mapping reads to 1st genome (done by JGI)
echo 'done'


# 2. filter umapped or poor mapped reads (none M < 30bp)
#samtools view $qry\.bam | grep -v ^@ | awk '$6=="151M" {print "@"$1"\n"$10"\n+\n"$11}'
python bam_flt.py foo.bam > foo.bam.flt.bam
bedtools bamtofastq -i foo.bam.flt.bam -fq /dev/stdout > foo.bam.flt.bam.fq


# 3. extract mapped range from 1st genome (bedtools)
bedtools merge -d 50 -i foo.bed > foo.bed.merge.bed
bedtools getfasta -fi genome.fasta -bed foo.bed.merge.bed > bed.fasta


# 4. using filtered reads to assembly (idba, megahit, platanus)
#platanus_allee assemble -f foo.bam.flt.bam.fq -o output
idba -r foo.bam.flt.bam.fq.fa -o output --pre-correction --mink 71 --maxk 101


# 5. using assembly from filtered reads to do gap closing (gatb, metaassembler)
# using metaassembler
./Metassemble_script.sh


# 6. construct dbg by cdbg_search
../a3.x inputfile outputfile kFile


# 7. identify FRs by FindFRs2
java -jar FindFRs.jar -d input.dot -f genome -a .3 -k 10



# snp and plink
# 1. merge all the vcf (vcf-tools)
# 2. convert the vcf to file required by plink (plink)
# 3. using plink and R to get qqnorm and manh plot (plink and R)
# 4. using gemma
gemma -bfile 2000 -gk 2 -o kin


# gene prediction
# 1. augustus model traning:
# (1). bam file
# (2). EST + genomethreader
# (3). gff file
# 2. de novo:
# (1). augustuts + trained model
# (2). genemark-es or et
# 3. homology prediction
# (1). exonerate + swissprot
# (2). genomethreader + swissprot
# 4. Utr prediction:
# (1) PASA
