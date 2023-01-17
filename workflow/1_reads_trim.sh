#!/bin/bash
## B73 reads trimming ##
java -Xmx16g -jar trimmomatic-0.38.jar PE \
./input/B73.R1.fq.gz ./input/B73.R2.fq.gz \
./cache/B73_trim.R1.fq.gz ./cache/B73_unpaired.R1.fq.gz \
./cache/B73_trim.R2.fq.gz ./cache/B73_unpaired.R2.fq.gz \
ILLUMINACLIP:./input/TruSeq3-PE.fa:3:20:10:1:true LEADING:3 TRAILING:3 SLIDINGWINDOW:4:13 MINLEN:40 2> ./cache/B73_trim.log

## A188 reads trimming ##
java -Xmx16g -jar trimmomatic-0.38.jar PE \
./input/A188.R1.fq.gz ./input/A188.R2.fq.gz \
./cache/A188_trim.R1.fq.gz ./cache/A188_unpaired.R1.fq.gz \
./cache/A188_trim.R2.fq.gz ./cache/A188_unpaired.R2.fq.gz \
ILLUMINACLIP:./input/TruSeq3-PE.fa:3:20:10:1:true LEADING:3 TRAILING:3 SLIDINGWINDOW:4:13 MINLEN:40 2> ./cache/A188_trim.log