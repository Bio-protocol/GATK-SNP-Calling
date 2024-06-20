#!/bin/bash
## gatk index for reference genome ##
samtool faidx ./input/B73Ref4.fa
gatk CreateSequenceDictionary -R ./input/B73Ref4.fa -O ./input/B73Ref4.dict

## GATK SNP Calling (2 samples together) ##
gatk HaplotypeCaller --java-options '-Xmx24g' -R ./input/B73Ref4.fa -I ./cache/B73.RG.RD.bam -I ./cache/A188.RG.RD.bam -O ./output/B73_A188.raw.0.vcf

## GATK Short Variant Calling (with bam list) ##
gatk HaplotypeCaller --java-options '-Xmx24g' \
-R ./input/B73Ref4.fa \
-I ./input/bam.list \
-O ./output/B73_A188.raw.0.vcf

