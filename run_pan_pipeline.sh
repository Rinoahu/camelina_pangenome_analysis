#!/bin/bash

#awk '{if($2<1000) print $1}' pfq.txt > allow.txt
#python pan_genome.py -l .006 -u .9 -i header.fsa -g all_pep.fsa.xyz.mcl.flt -r allow.txt &>header.fsa.pan

python pan_genome.py -l .01 -u .9 -i header.fsa -g all_pep.fsa.xyz.mcl.flt &>header.fsa.pan

grep \^test header.fsa.pan > tmp.txt

Rscript R_script.rs
