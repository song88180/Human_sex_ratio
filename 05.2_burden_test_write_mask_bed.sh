#!/bin/bash
# regenie_masks.txt contains only one line "M3 LoF,missense(5/5)"
# --phenoFile and --phenoCol does not matter, just to ensure that the script can run without error.
# ukb23158_500k_OQFE.90pct10dp_qc_variants.txt,  ukb23158_500k_OQFE.annotations.txt.gz, 
# and ukb23158_500k_OQFE.sets.txt.gz are provided by UKB Research Analysis Platform.

N=$1

regenie \
  --step 2 \
  --bed ./ukb23158_c${N}_b0_v1 \
  --ignore-pred \
  --skip-test \
  --exclude ./ukb23158_500k_OQFE.90pct10dp_qc_variants.txt \
  --phenoFile ./GWAS_adj_extended_new.phe \
  --phenoCol N_full_brothers \
  --anno-file ./ukb23158_500k_OQFE.annotations.txt.gz \
  --set-list ./ukb23158_500k_OQFE.sets.txt.gz \
  --mask-def ./regenie_masks.txt \
  --aaf-bins 0.01 \
  --bsize 300 \
  --threads 16 \
  --write-mask \
  --nauto 24 \
  --out ./chr${N}_burden
