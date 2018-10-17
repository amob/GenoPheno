#!/bin/bash -l

module load gcc jdk/1.8 tassel/5 vcftools/0.1.13

run_pipeline.pl -Xmx15g -fork1 -h5 AOMexGBS2_ZeaGBSv27impV5.h5 -excludeTaxa BLANK:250481877,amA:250481966,amB:250481967,amC:250481968,amD:250481969,amE:250481970 -export AOMexGBS2_ZeaGBSv27impV5 -exportType VCF -runfork1

vcftools --vcf AOMexGBS2_ZeaGBSv27impV5.vcf --out missing_AOMexGBS2 --missing-indv

vcftools --vcf MinimallyFilteredWithNAM_InbredLandraces_InbredTeosintes.vcf --out missing_swarts2017  --missing-indv

srun R CMD BATCH 09_keeptaxa.R
#THEN SEE R FILE FOR VARIOUS STEPS REQUIRED TO CONVERT COORDINATES OF SNPS TO SAME REFERENCE MAP BEFORE NEXT LINE







