

swarts <- read.table("~/GenoPheno/MinimallyFilteredWithNAM_InbredLandraces_InbredTeosintes.vcf",header=F,stringsAsFactors=F)
swartspos <- paste(swarts[,1],swarts[,2],sep=".")
annamex <- read.table("~/GenoPheno/AOMexGBS2_ZeaGBSv27impV5.vcf",header=F,stringsAsFactors=F)
#swarts snps are in refgenv3 (annamex are v2)
#automatically skips all # out rows, so will have to add in the colnames line later.
annapos <- paste(annamex[,1],annamex[,2],sep=".")
rm(annamex)

options(scipen = 999)

swartsV3toV2 <- read.table("~/GenoPheno/output_MinimallyFilteredWithNAM_InbredLandraces_InbredTeosintes_forconversion.bed",header=F,stringsAsFactors=F)
swartsV2toV3 <- read.table("~/GenoPheno/output_output_MinimallyFilteredWithNAM_InbredLandraces_InbredTeosintes_forconversion.bed",header=F,stringsAsFactors=F)

wasconv.pos <- paste(swartsV2toV3[,1],swartsV2toV3[,2],sep=".")
conv.swarts <- swarts[swartspos%in%wasconv.pos,]
rm(swarts)

conv.swarts[,1:2] <- swartsV3toV2[,1:2]
conv.swarts[,3] <- paste("S",conv.swarts[,1],"_",conv.swarts[,2],sep="")
conv.swarts.order <- conv.swarts[order(conv.swarts[,1],conv.swarts[,2]),]
write.table(conv.swarts.order,"~/GenoPheno/MinimallyFilteredWithNAM_InbredLandraces_InbredTeosintes_agpv3_ConvertedTo_agpv2.vcf",row.names=F,quote=F,col.names=F,sep="\t")


