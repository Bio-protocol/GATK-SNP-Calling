![image](https://github.com/Bio-protocol/GATK-SNP-Calling/assets/54077130/a50034fe-3804-45d3-9c8c-17df6a1f8fa1)# GATK Variant Discovery Pipeline
## Introduction
Phenotypic variations of most biological traits are largely driven by genomic variants. The single nucleotide polymorphism (SNP) is the most common form of genomic variants. Multiple algorithms have been developed for discovering variants, including SNPs, with next generation sequencing (NGS) data. Here we present a widely used variant discovery pipeline based on the software Genome Analysis ToolKits (GATK). The pipeline uses whole genome sequencing (WGS) data as input data and includes read mapping, variant calling, and variant filtering processes. This pipeline has been successfully applied to many genomic projects and represent a solution for variant calling using NGS data.
 
Here we create a protocol to do SNP Calling for Whole Genome Sequencing (WGS) data, including <br />
(1) Align reads to reference genome with [BWA](http://bio-bwa.sourceforge.net/) <br />
(2) SNP Calling based on [GATK4 Best Practices Workflows](https://gatk.broadinstitute.org/hc/en-us/articles/360035535932-Germline-short-variant-discovery-SNPs-Indels-) <br />
(3) Provide the suggesting filtering parameters for SNPs filtration <br />

1. **graphs**: The graphs/figures produced during the analysis.
2. **input**: The raw input data.
3. **lib**: The source code, functions, or algorithms used within the workflow.
4. **cache**: The intermediate files when running the workflow.
5. **output**: The final output results of the workflow.
6. **workflow**: Step by step pipeline.
## Overview of the workflow
This is the workflow to show a step-by-step pipeline to call SNPs.
![](https://github.com/hecheng90/GATK-SNP-Calling/blob/main/graphs/Diagram.png)
## Installation
### Required software and installation

All the softwares were installed on Linux (x64) operating system (OS).

#### Installing [Anaconda](https://www.anaconda.com/download)
```
wget https://repo.anaconda.com/archive/Anaconda3-2024.02-1-Linux-x86_64.sh
bash Anaconda3-2024.02-1-Linux-x86_64.sh
echo 'export PATH="~/anaconda3/bin:$PATH"' >> ~/.bashrc
source ~/.bashrc
```

#### All the softwares can be installed through the GATK_SNP.yaml file in the github repository.
> One-command run to install all softwares by GATK_SNP.yaml 
```
conda env create -f GATK_SNP.yaml
source activate GATK_SNP
```

#### (Optional) You can also install softwares one by one manually
> Create and load conda environment
```
conda create -n GATK_SNP
source activate GATK_SNP
```
#### Installing [Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic)
> Software for reads trimming
```
conda install -c bioconda trimmomatic
```
#### Installing [BWA](http://bio-bwa.sourceforge.net/)
> Software for reads alignment
```
conda install -c bioconda bwa
```
#### Installing [SAMtools](http://www.htslib.org/)
> Software for SAM/BAM file treatment
```
conda install -c bioconda samtools
```
#### Installing [GATK4](https://gatk.broadinstitute.org/hc/en-us)
> Software for SNP Calling
```
conda install -c bioconda gatk4
```

#### Installing [vcftools](http://vcftools.sourceforge.net/)
> Software to measure vcf files
```
conda install -c bioconda vcftools
```

## Input Data
The sample data for GATK SNP Calling workflow include: 
- B73/A188.R1/2.fq.gz: Paired-end whole genome sequencing reads for two test samples.
- B73Ref4.fa: The maize B73 version4 reference genome.
- TruSeq3-PE.fa: The adaptor sequences for reads trimming
The raw reads and adaptor sequences were included in input folder. The B73 version4 reference genome can be downloaded from MaizeGDB (https://download.maizegdb.org/Zm-B73-REFERENCE-GRAMENE-4.0/ Zm-B73-REFERENCE-GRAMENE-4.0.fa.gz)

## Major Step
### Step0: Download necessary files
> You can download the github repository with this command line, which includes all the files and scripts need for this protocol:
```
git clone https://github.com/Bio-protocol/GATK-SNP-Calling.git
cd GATK-SNP-Calling
```

### Step1: Raw reads trimming
> Note: The trimmomatic software should run on 1.8.0+ Java environment
```
sh ./workflow/1_reads_trim.sh
```
- Script content of 1_reads_trim.sh
```
#!/bin/bash
## B73 reads trimming ##
trimmomatic PE \
./input/B73.R1.fq.gz ./input/B73.R2.fq.gz \
./cache/B73_trim.R1.fq.gz ./cache/B73_unpaired.R1.fq.gz \
./cache/B73_trim.R2.fq.gz ./cache/B73_unpaired.R2.fq.gz \
ILLUMINACLIP:./input/TruSeq3-PE.fa:3:20:10:1:true LEADING:3 TRAILING:3 SLIDINGWINDOW:4:13 MINLEN:40 2> ./cache/B73_trim.log

## A188 reads trimming ##
trimmomatic PE \
./input/A188.R1.fq.gz ./input/A188.R2.fq.gz \
./cache/A188_trim.R1.fq.gz ./cache/A188_unpaired.R1.fq.gz \
./cache/A188_trim.R2.fq.gz ./cache/A188_unpaired.R2.fq.gz \
ILLUMINACLIP:./input/TruSeq3-PE.fa:3:20:10:1:true LEADING:3 TRAILING:3 SLIDINGWINDOW:4:13 MINLEN:40 2> ./cache/A188_trim.log
```
## Step2: BWA reads alignment
```
sh ./workflow/2_bwa_align.sh
```
- Script content of 2_bwa_align.sh
```
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
```
### Parameters for BWA:
> - **-t INT**: Number of threads
> - **-R STR**: Complete read group header line. ’\t’ can be used in STR and will be converted to a TAB in the output SAM. The read group ID will be attached to every read in the output. An example is ’@RG\tID:foo\tSM:bar’.
## Step3: SAM filtering and compressing
> The parameters used here are set for test sample alignment
```
sh ./workflow/3_SAM_filter.sh
```
- Script content of 3_SAM_filter.sh
```
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
```
### Parameters for samparser.bwa.pl:
> - **--input|i:** SAM file
> - **--identical|e:** minimum matched and identical base length, default=30bp
> - **--mismatches|mm|m:** two integers to specify the number of mismatches out of the number of basepairs of the matched region of reads; (matched regions are not identical regions, mismatch and indel could occur e.g., --mm 2 36 represents that <=2 mismatches out of 36 bp
> - **--tail:** the maximum bp allowed at each side, two integers to specify the number of tails out of the number of basepairs of the reads, not including "N", e.g., --tail 3 75 represents that <=3 bp tails of 75 bp of reads without "N"
> - **--gap:** if a read is split, the internal gap (bp) allowed, default=5000bp
> - **--mappingscore:** the minimum mapping score, default=40;
> - **--insert:** insert range, e.g., 100 600 (default)
> - **--help:** help information

## Step4: Data Pre-processing for GATK SNP Calling
> GATK and Java environment are required for this step
```
sh ./workflow/4_data_preprocess.sh
```
- Script content of 4_data_preprocess.sh
```
#!/bin/bash
## (Optinal) Add labels for inputs ##
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
```
> **(Optional) Recalibrate Base Quality Scores**
> ```
> You can use known variation resources to improve SNP Calling
> ## Generates recalibration table for Base Quality Score Recalibration (BQSR) ##
> gatk BaseRecalibrator \
> -I input.RG.RD.bam \
> -R reference.fasta \
> --known-sites sites_of_variation.vcf \
> --known-sites another/optional/setOfSitesToMask.vcf \
> -O recal_data.table
>
> ## Apply base quality score recalibration ##
> gatk ApplyBQSR \
> -R reference.fasta \
> -I input.RG.RD.bam \
> --bqsr-recal-file recalibration.table \
> -O output.RG.RD.BQSR.bam    
> ```
## Step5: GATK SNP Calling
> Java environment (1.8.0+) is required for this step
```
sh ./workflow/5_GATK.sh
```
- Script content of 5_GATK.sh
```
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
```
### Notes:
> (1)Here we use the default parameters of GATK HapotypeCaller to do SNP Calling. The details of each GATK parameters can be found [here](https://gatk.broadinstitute.org/hc/en-us/articles/360036853331-HaplotypeCaller).
> (2)We suggest to do GATK Short Variant Calling with all samples together. If there are too many samples to load at the same time, you can use a bam list file including all sample information instead (## GATK Short Variant Calling (with bam list) ##)

## Step6: SNPs Filtering
```
sh ./workflow/6_SNP_filter.sh
```
- Script content of 6_SNP_filter.sh
```
#!/bin/bash
## select bi-allelic SNPs ##
vcf=./output/B73_A188.raw.0.vcf
ref=./input/B73Ref4.fasta
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
-O ./output/B73_A188.HF.2.vcf.gz

## Extract PASS SNPs ##
gatk SelectVariants --java-options '-Xmx24g' -R $ref -V ./output/B73_A188.HF.2.vcf.gz \
--exclude-filtered \
-O ./output/B73_A188.PASS.3.vcf.gz

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
```
> The recommended hard-filtering parameters for SNPs by GATK can be found [here](https://gatk.broadinstitute.org/hc/en-us/articles/360035890471-Hard-filtering-germline-short-variants).

## Expect Results
### Reads Trimming
The reads trimming information were saved in ./cache/*_trim.log files.
```
Input Read Pairs: 200000 
Both Surviving: 200000 (100.00%) 
Forward Only Surviving: 0 (0.00%)
Reverse Only Surviving: 0 (0.00%) Dropped: 0 (0.00%)
```
### Reads Alignment and Filtering
The reads alignment and filtering information were saved in ./cache/*.parse.log files.
```
./cache/B73.sam Total reads in the SAM output   200000
./cache/B73.sam Reads could be mapped   199929
./cache/B73.sam Passing criteria reads  140931
./cache/B73.sam Unmapped reads  969
```
### Remove duplicates
The bam remove duplicates information were saved in ./cache/*.metrics.txt files.
| LIBRARY | UNPAIRED_READS_EXAMINED | READ_PAIRS_EXAMINED | SECONDARY_OR_SUPPLEMENTARY_RDS | UNMAPPED_READS | UNPAIRED_READ_DUPLICATES | READ_PAIR_DUPLICATES | READ_PAIR_OPTICAL_DUPLICATES | PERCENT_DUPLICATION | ESTIMATED_LIBRARY_SIZE |
| ------- | ----------------------- | ------------------- | ------------------------------ | -------------- | ------------------------ | -------------------- | ---------------------------- | ------------------- | ---------------------- |
| A | 0 | 275194 | 0 | 0 | 0 | 173 | 37 | 0.000629 | 278258915 |
### GATK SNP Calling
The GATK SNP Calling results were saved in ./output/*.vcf files.
| CHROM | POS | ID | REF | ALT | QUAL | FILTER | INFO | FORMAT | B73 | A188 |
| ----- | --- | -- | --- | --- | ---- | ------ | ---- | ------ | ------- | ------- |
| 1 | 242449 | . | T | C | 23.19 | PASS | AC=2;AF=1.00;AN=2;DP=2;ExcessHet=3.0103;FS=0.000;MLEAC=2;MLEAF=1.00;MQ=46.50;QD=11.60;SOR=0.693 GT:AD:DP:GQ:PL | GT:AD:DP:GQ:PL | 1/1:0,2:2:6:49,6,0 | ./. |
> The VCF format information can be found in [VCF format wiki](https://en.wikipedia.org/wiki/Variant_Call_Format)
## License
It is a free and open source software, licensed under [GPLv3](https://github.com/github/choosealicense.com/blob/gh-pages/_licenses/gpl-3.0.txt)
