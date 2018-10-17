#!/bin/bash -l

module load gcc jdk/1.8 tassel/5.2.14 plink/1.90 vcftools/0.1.13

#fixing header
head -n 835 MinimallyFilteredWithNAM_InbredLandraces_InbredTeosintes.vcf > header.txt
cat header.txt MinimallyFilteredWithNAM_InbredLandraces_InbredTeosintes_agpv3_ConvertedTo_agpv2.vcf >  MinimallyFilteredWithNAM_InbredLandraces_InbredTeosintes_agpv3_ConvertedTo_agpv2_newheader.vcf

run_pipeline.pl -Xmx15g -fork1 -vcf AOMexGBS2_ZeaGBSv27impV5.vcf -includeTaxaInFile AOMexGBS2_indFilt70.txt -fork2 -vcf MinimallyFilteredWithNAM_InbredLandraces_InbredTeosintes_agpv3_ConvertedTo_agpv2_newheader.vcf -includeTaxaInFile SwartsKeepTaxa_indFilt70.txt  -combine3 -input1 -input2 -mergeGenotypeTables -fork4 -filterAlign -input3 -filterAlignMinCount 421 -export MergeSwartsAnna_80 -exportType Plink -runfork1 -runfork2 -runfork3 -runfork4 
#421 is 80% of 526 samples. sites with 20% or less missing data

run_pipeline.pl -Xmx15g -fork1 -Plink -ped MergeSwartsAnna_80.plk.ped -map MergeSwartsAnna_80.plk.map -PrincipalComponentsPlugin -covariance true -endPlugin -export pca_MergeSwartsAnna_80.txt -runfork1







