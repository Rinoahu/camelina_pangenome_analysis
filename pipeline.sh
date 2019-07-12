#!/bin/bash

########################################################################################
# the first sequenced genome of camelina
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
