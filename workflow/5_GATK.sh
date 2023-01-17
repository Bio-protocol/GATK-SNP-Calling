#!/bin/bash
## gatk index for reference genome ##
samtool faidx ./input/B73Ref4.fasta
mv B73Ref4.fasta.fai ./input/
gatk CreateSequenceDictionary -R ./input/B73Ref4.fasta -O ./input/B73Ref4.dict

## GATK SNP Calling (2 samples together) ##
gatk HaplotypeCaller --java-options '-Xmx24g' -R ./input/B73Ref4.fasta -I ./cache/B73.RG.RD.bam -I ./cache/A188.RG.RD.bam -O ./output/B73_A188.raw.0.vcf

## GATK Short Variant Calling (with bam list) ##
gatk HaplotypeCaller --java-options '-Xmx24g' \
-R ./input/B73Ref4.fasta \
-I ./input/bam_list.txt \
-O ./output/B73_A188.raw.0.vcf

