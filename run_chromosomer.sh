#!/bin/bash

# the reference
ref=$1

# the contigs
qry=$2

echo "get length of fasta"
chromosomer fastalength $qry $qry\.length


echo "blastn alignment"
makeblastdb -in $ref -dbtype nucl


blastn -query $qry -db $ref -outfmt 6 -out $qry\_alignments.txt


chromosomer fragmentmap $qry\_alignments.txt 100 $qry\.length $qry\_map.txt

echo "assembly"
chromosomer assemble $qry\_map.txt $qry assembled_chromosomes.fa





