#!/bin/sh
# an shell for orthomcl pipeline

all_in_1=$1
# access the mysql database
#mysql -u root -pv2110<<EOF
#drop database orthomcl;
#create database orthomcl;
#EOF

#mysql -u root -pv2110<<EOF
#create database orthomcl;
#EOF


#cp /home/xh/gnu_tools/orthomclSoftware-v2.0.9/orthomcl.config ./
rm -rf orthomcl_sqliteTEST*
orthomclInstallSchema orthomcl.config
#rm -rf compliantFasta
#mkdir compliantFasta
#cp $1 compliantFasta/
#cd compliantFasta
#orthomclAdjustFasta fungi ../$all_in_1 fungi
#cd ..
#orthomclFilterFasta compliantFasta 10 20

#exit 0

# use blast or blat to do self-blast search
#formatdb -i goodProteins.fasta -p T
#/usr/bin/time -v legacy_blast blastall -p blastp -i goodProteins.fasta -o goodProteins.fasta.blast -d goodProteins.fasta -m 8 -a 16 -F T -e 1e-5 -v 1000000 -b 1000000

#diamond=/home/xh/Download/tools/diamond/diamond
#$diamond makedb --in goodProteins.fasta -d goodProteins.fasta
#/usr/bin/time -v $diamond blastp --quiet -q goodProteins.fasta -d goodProteins.fasta --more-sensitive -k 1000000 -e 1e-5 -o goodProteins.fasta.m8

#exit 1

#ln -sf goodProteins.fasta.blast8 my_blast_results
ln -sf goodProtein_new.sc my_blast_results
#ln -s /home/hx/tools/new_genome_project/antismash2.0/antismash-2.0.1/kaks/../blastout/goodProteins.fasta.blast my_blast_results
rm -rf pairs

orthomclBlastParser my_blast_results compliantFasta > similarSequences.txt
perl -p -i -e 's/\t(\w+)(\|.*)orthomcl/\t$1$2$1/' similarSequences.txt
perl -p -i -e 's/0\t0/1\t-357/' similarSequences.txt

orthomclLoadBlast orthomcl.config similarSequences.txt
orthomclPairs orthomcl.config log_file cleanup=yes
orthomclDumpPairsFiles orthomcl.config
mcl mclInput --abc -I 1.5 -o mclOutput -te 8
cat mclOutput | orthomclMclToGroups group 00001 > groups.txt

