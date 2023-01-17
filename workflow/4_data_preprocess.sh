#!/bin/bash
## Add labels for inputs ##
gatk AddOrReplaceReadGroups --java-options '-Xmx24g' -I ./cache/B73.parse.sort.bam -O ./cache/B73.RG.bam -LB A -PL illumina -PU A -SM B73 -ID B73
gatk AddOrReplaceReadGroups --java-options '-Xmx24g' -I ./cache/A188.parse.sort.bam -O ./cache/A188.RG.bam -LB A -PL illumina -PU A -SM A188 -ID A188
#parameters for AddOrReplaceReadGroups#
#-LB: Read-Group library#
#-PL: Read-Group platform#
#-PU: Read-Group platform unit#
#-SM: Read-Group sample name#
#-ID: Read-Group ID#

## MarkDuplicates ##
gatk MarkDuplicates --java-options '-Xmx24g' --REMOVE_DUPLICATES true -I ./cache/B73.RG.bam -M ./cache/B73.metrics.txt -O ./cache/B73.RG.RD.bam
gatk MarkDuplicates --java-options '-Xmx24g' --REMOVE_DUPLICATES true -I ./cache/A188.RG.bam -M ./cache/A188.metrics.txt -O ./cache/A188.RG.RD.bam

## BAM index ##
cd cache
samtools index B73.RG.RD.bam
samtools index A188.RG.RD.bam
rm *.parse.sort.bam
rm *.RG.bam