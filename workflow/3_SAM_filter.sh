#!/bin/bash
## SAM filtering ##
perl ./lib/samparser.bwa.pl -i ./cache/B73.sam -e 100 -m 3 100 --tail 5 100 --gap 0 --insert 100 800 1>./cache/B73.parse.sam 2>./cache/B73.parse.log
perl ./lib/samparser.bwa.pl -i ./cache/A188.sam -e 100 -m 3 100 --tail 5 100 --gap 0 --insert 100 800 1>./cache/A188.parse.sam 2>./cache/A188.parse.log

## BAM sort & index ##
cd cache
samtools view -u B73.parse.sam | samtools sort -o B73.parse.sort.bam
samtools view -u A188.parse.sam | samtools sort -o A188.parse.sort.bam
samtools index B73.parse.sort.bam
samtools index A188.parse.sort.bam
rm *.sam
cd ../
