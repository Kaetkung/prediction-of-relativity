#!/bin/bash

# Convert VCF to PLINK binary format
./plink --vcf start.vcf.gz --make-bed --vcf-half-call missing --out df1

# List duplicate variants
./plink --bfile df1 --list-duplicate-vars ids-only suppress-first

# Perform quality control and generate compressed VCF
./plink --bfile df1 --geno 0.05 --mind 0.05 --maf 0.05 --autosome --snps-only just-acgt --exclude plink.dupvar --allow-no-sex --nonfounders --recode vcf-iid bgz --out afterQC

# Index the compressed VCF file
bcftools index afterQC.vcf.gz

# Split VCF file by chromosome
for i in {1..22}; do
  bcftools view -O z -o input_chr${i}.vcf.gz -r ${i} afterQC.vcf.gz
done

# Phase each chromosome using Beagle without imputation
for i in {1..22}; do
  java -Xmx12g -Xms12g -jar beagle.22Jul22.46e.jar gt=input_chr${i}.vcf.gz ref=chr${i}.1kg.phase3.v5a.b37.bref3 map=plink.chr${i}.GRCh37.map impute=false out=phased_chr${i}_output
done

# Index phased VCF files
for i in {1..22}; do
  bcftools index phased_chr${i}_output.vcf.gz
done

# Concatenate phased chromosome files
bcftools concat -O z -o phased_output_combined.vcf.gz phased_chr{1..22}_output.vcf.gz

# Convert VCF to PLINK binary format
./plink --vcf phased_output_combined.vcf.gz --make-bed --vcf-half-call missing --out df2

# Estimate pairwise IBD
./plink --bfile df2 --genome --min 0.015 --memory 12582912

# concat genetic map
cat plink.chr*.GRCh37.map > combined_genetic_map.map

# generate segment
java -Xmx12g -Xms12g -jar refined-ibd.17Jan20.102.jar gt=phased_output_combined.vcf.gz map=combined_genetic_map.map out=ibd_segments length=7

# merge IBD
gunzip -c ibd_segments.ibd.gz | java -Xmx12g -Xms12g -jar merge-ibd-segments.17Jan20.102.jar phased_output_combined.vcf.gz combined_genetic_map.map 0.6 1 > merged_ibd_segments.ibd
