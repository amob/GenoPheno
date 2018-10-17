#!/bin/bash -l

module load gcc jdk/1.8 tassel/5

run_pipeline.pl -Xmx3g -fork1 -h5 AOMexGBS2_ZeaGBSv27impV5.h5 -excludeTaxa BLANK:250481877,amA:250481966,amB:250481967,amC:250481968,amD:250481969,amE:250481970 -fork2 -filterAlign -input1 -filterAlignMinCount 86 -export AOMexGBS2_ZeaGBSv27impV5_95filter -exportType Plink -runfork1 -runfork2

#filtering sites for 86 taxa with info out of 90 taxa, 86 out of 90 is 95%
#The taxa in the exclude taxa list are a population not included in the experiment, and the blank.

