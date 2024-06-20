#!/bin/bash
## select bi-allelic SNPs ##
vcf=./output/B73_A188.raw.0.vcf
ref=./input/B73Ref4.fa
gatk SelectVariants \
        -R $ref \
        -V $vcf \
        --restrict-alleles-to BIALLELIC \
        -select-type SNP \
        -O ./output/B73_A188.bi.1.vcf.gz &>./output/B73_A188.bi.log

## Hard-filtering germline short variants ##
gatk VariantFiltration --java-options '-Xmx24g' -R $ref -V ./output/B73_A188.bi.1.vcf.gz \
--filter-expression "QD < 2.0 || MQ < 40.0 || FS > 60.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0 || SOR > 3.0" \
--filter-name "hard_filter"  \
-O ./output/B73_A188.HF.2.vcf.gz &>./output/B73_A188.HF.log

## Extract PASS SNPs ##
gatk SelectVariants --java-options '-Xmx24g' -R $ref -V ./output/B73_A188.HF.2.vcf.gz \
--exclude-filtered \
-O ./output/B73_A188.PASS.3.vcf.gz &>./output/B73_A188.PASS.log

## (Optional) Convert heterozygous SNPs to Missing ##
gatk VariantFiltration --java-options '-Xmx24g' \
-V ./output/B73_A188.PASS.3.vcf.gz \
-O ./output/B73_A188.mark.hetero.4.vcf.gz \
--genotype-filter-expression "isHet == 1" \
--genotype-filter-name "isHetFilter"

gatk SelectVariants --java-options '-Xmx24g' SelectVariants \
-V ./output/B73_A188.mark.hetero.4.vcf.gz \
--set-filtered-gt-to-nocall \
-O ./output/B73_A188.het2miss.4.vcf.gz

## (Optional) Filter SNPs by Minor Allele Frequency (MAF) & Missing Rate (MR) ##
vcftools ./output/B73_A188.het2miss.4.vcf.gz \
--maf 0.02 --max-missing 0.2 --recode --recode-INFO-all \
--out ./output/B73_A188.MAF.MR.5.vcf.gz
