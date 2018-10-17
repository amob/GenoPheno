#install.packages("~/tarprograms/RAFM_1.2.tar.gz",repos=NULL) #these must be downloaded and installed, local directories will differ
#install.packages("~/tarprograms/SparseM_1.7.tar.gz",repos=NULL) 
#install.packages("~/tarprograms/corpcor_1.6.8.tar.gz",repos=NULL)

library(RAFM)
##RAFM requires genotype data in the format where the first column is the population, columns 2&3 have alleles 1 and 2 for locus 1; 4&5 alleles 1 and 2 for locus 2.
genodata <- read.table("~/GenoPheno/AOMexGBS2_ZeaGBSv27impV5_95filter.plk.ped")#first 6 columns of this file are not locus data. 
genodata <- genodata[,-c(1,3:6)]#. rm columns that are not locus data

####SUBSET TO A RANDOM SET
poss <- seq(from=2,to=(ncol(genodata)-1),by=2) #for all loci in dataset, even numbers representing allele 1.
##choose random subset of poss
tenK <-sort(sample(poss,10000))#the 10k loci
colnums<-c(1,unlist(lapply(tenK, function(z) c(z,z+1))))#the column numbers needed
genosub <- genodata[,colnums]

pmiss2 <- sapply(1:nrow(genodata), function(z) length(which(genodata[z,2:ncol(genodata)]==0)) ) /length(2:ncol(genodata))
mean(pmiss2[-c(25,75)])#removing taxa that will not be used, average missing data of an individaul plant genotype data: 0.0088
pmiss <- sapply(1:nrow(genosub), function(z) length(which(genosub[z,]==0)))#missing data of individuals in genosub dataset
genosub <- genosub[-(which((pmiss/20000)>.5)),]#two columns per locus: allele 1 and allele 2, divide by 20k not 10k
#this removes two samples with very high missing data
rm(genodata)

genoid <- sapply(1:nrow(genosub), function(z) strsplit(as.character(genosub[z,1]),":")[[1]][1])
genoidtopop <- read.csv("~/GenoPheno/popmomlabel.csv",header=T,stringsAsFactors=F)
genopop <- sapply(1:length(genoid), function(z) genoidtopop$pop[which(genoidtopop$label==genoid[z])])
genopop[which(genopop=="M")] <-"ML"#"ml"#fixing malinalco so that has 2 letter abbreviation like all other populations
genosub[,1] <- genopop#write over first column of lane #s with genopop.
write.table(genosub,"~/GenoPheno/AOMexGBS2_ZeaGBSv27impV5_95filter_Rthin.ped",row.names=F,quote=F)



####RUN coancestry analysis
genosub <- read.table("~/GenoPheno/AOMexGBS2_ZeaGBSv27impV5_95filter_Rthin.ped",header=T)
all.coan <- do.all(genosub,20000,10000,10)#these are recommended mcmc inputs.
save(all.coan,file="allcoan_Rthin_f95_LG.txt")

####FIGURE
load("~/GenoPheno/allcoan_Rthin_f95_LG.txt")
torder <- c(10,5,9,7,3,4,6,1,8,2) #in the file as provided, they are in order that were given in the original file (i.e. alphabetical, except that malinalco is before matias cuijingo)
envmatall <-read.csv("~/GenoPheno/AlphabeticalPopEnvDat.csv")
templabels <- as.character(envmatall$TAnn[torder]/10) # tz, ml, tx, mt, da, fp, mc, cl, tc, cu
templabels[4] <- "15.0"
templabels[9] <- "13.0"

library(gplots)
library(fields)
library(SDMTools)

thetaAll <- apply(all.coan$theta,c(1,2),mean)
colnames(thetaAll) <- c("CL","CU","DA","FP","M","MC","MT","TC","TX","TZ")
thetaAlln <- thetaAll[torder,torder]
thetaAlln[upper.tri(thetaAlln,diag=F)] <- NA

bkw<-colorRampPalette(c(rgb(1,1,1),rgb(0,0,0)))

pdf("~/GenoPheno/ALLrafmcoan_f95_LG2MXANNA.pdf",width=4,height=4)
image.plot(seq(from=.5,to=9.5,by=1),seq(from=.5,to=9.5,by=1),t(thetaAlln),axes=FALSE, legend.shrink=.7,
 ylab = "Mean Annual Temperature", xlab = "",cex.lab=1, ylim=c(0,10),col=bkw(10))
# the following pointless command is necessary to make the custom axis labels non-transparent, google revealed this among a number of other workarounds.
points(0,0)# now these will display properly
axis(side=2, at=seq(from=.5,to=9.5,by=1), labels=templabels, las=2, cex.axis = 0.8)
axis(side=3, at=seq(from=.5,to=9.5,by=1), labels=templabels, las=2, cex.axis = 0.8)
mtext("Coancestry", side=3,line=-1,adj=2)
dev.off()

