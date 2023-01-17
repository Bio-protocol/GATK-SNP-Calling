#!/bin/bash
## BWA index of reference genome ##
cd input
wget https://download.maizegdb.org/Zm-B73-REFERENCE-GRAMENE-4.0/Zm-B73-REFERENCE-GRAMENE-4.0.fa.gz
gunzip Zm-B73-REFERENCE-GRAMENE-4.0.fa.gz
mv Zm-B73-REFERENCE-GRAMENE-4.0.fa B73Ref4.fa
bwa index B73Ref4.fa

## BWA alignment ##
bwadb=./input/B73Ref4.fa
B73_R1=./cache/B73_trim.R1.fq.gz
B73_R2=./cache/B73_trim.R2.fq.gz
A188_R1=./cache/A188_trim.R1.fq.gz
A188_R2=./cache/A188_trim.R2.fq.gz
bwa mem -t 8 -R '@RG\tID:B73\tSM:B73' $bwadb $B73_R1 $B73_R2 > ./cache/B73.sam
bwa mem -t 8 -R '@RG\tID:A188\tSM:A188' $bwadb $A188_R1 $A188_R2 > ./cache/A188.sam