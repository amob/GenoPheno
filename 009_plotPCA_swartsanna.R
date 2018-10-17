#r script to make a plot of pca of a subset of swarts and  present teosinte populationss, to see if low elevation ones cluster with parviglumis

pcaSA <- read.table("~/GenoPheno/pca_MergeSwartsAnna_801.txt",header=T,skip=2,stringsAsFactors=F)
pcaSA.2 <- read.table("~/GenoPheno/pca_MergeSwartsAnna_802.txt",header=F,skip=1,stringsAsFactors=F)
pcaSA.3 <- read.table("~/GenoPheno/pca_MergeSwartsAnna_803.txt",header=T,stringsAsFactors=F)

SampleData <- read.csv("~/GenoPheno/Swarts_TableS3.csv",header=T,stringsAsFactors=F)

pcaSA$subspecies <- sapply(pcaSA[,1],function(z) ifelse(z%in%SampleData[,1],SampleData[which(SampleData[,1]==z),4],"presumedmex"))
#474 to 482 are currently malinalco
#518 to 526 are currently tepoztlan
pcaSA$subspecies[c(474:482,518:526)] <- "lowEmex"

colvector <- (as.factor(pcaSA$subspecies))
levels(colvector) <- c("green","black","forestgreen","blue","darkgreen")
#> levels(as.factor(pcaSA$subspecies))
#[1] "lowEmex"     "mays"        "mexicana"    "parviglumis" "presumedmex"
colvector <- as.character(colvector)

percentages <- paste(c("pca1 ","pca2 ","pca3 ","pca4 ","pca5 "),as.character(round(pcaSA.2[1:5,3]*100,digits=2)),rep("%",times=5),sep="")

pdf("~/GenoPheno/swartsannaPCA.pdf")
layout(matrix(1:4,ncol=2,byrow=T))
plot(pcaSA[,3]~pcaSA[,2],col=colvector,pch=1,ylab=percentages[2], xlab=percentages[1])
plot(pcaSA[,4]~pcaSA[,3],col=colvector,pch=1,ylab=percentages[3],xlab=percentages[2])
legend(-33,-10,c("mays","parviglumis","mexicana","high elevation","low elevation"),fill = c("black","blue","forestgreen","darkgreen","green"),bty="n")
plot(pcaSA[,5]~pcaSA[,4],col=colvector,pch=1,ylab=percentages[4],xlab=percentages[3])
plot(pcaSA[,6]~pcaSA[,5],col=colvector,pch=1,ylab=percentages[5],xlab=percentages[4])
dev.off()
