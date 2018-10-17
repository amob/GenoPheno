

swarts <- read.table("~/GenoPheno/MinimallyFilteredWithNAM_InbredLandraces_InbredTeosintes.vcf",header=F,stringsAsFactors=F)
swartspos <- paste(swarts[,1],swarts[,2],sep=".")

options(scipen = 999)
bedexportswarts <- cbind(swarts[,1],swarts[,2],(swarts[,2]+1))
write.table(bedexportswarts,"~/GenoPheno/MinimallyFilteredWithNAM_InbredLandraces_InbredTeosintes_forconversion.bed",quote=F,row.names=F,col.names=F)
