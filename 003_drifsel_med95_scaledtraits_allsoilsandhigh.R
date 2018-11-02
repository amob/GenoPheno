#install.packages("RAFM_1.2.tar.gz",repos=NULL)
#install.packages("SparseM_1.7.tar.gz",repos=NULL)
#install.packages("corpcor_1.6.8.tar.gz",repos=NULL)
#install.packages("driftsel_2.1.2.tar.gz",repos=NULL)
 

library(driftsel)
library(RAFM)
load("~/GenoPheno/allcoan_Rthin_f95_LG.txt")#object name will be all.coan

#read trait data
allpheno <- read.csv("~/GenoPheno/phenomat.csv",stringsAsFactors=F,header=T)
reds <- (allpheno$red)/(allpheno$red + allpheno$blue + allpheno$green)#colors as ratios with total brightness
greens <- (allpheno$green)/(allpheno$red + allpheno$blue + allpheno$green)
blues <- (allpheno$blue)/(allpheno$red + allpheno$blue + allpheno$green)

allpheno <- allpheno[,-c(3,10:12,14,16)]#cuts male and female flower day difference, number of leaves, flowering leaves, length leaves, raw redness, and raw blueness traits
allpheno$green <- greens #replace raw green measure with ratio of green

#pedigree file
allped <- read.csv("~/GenoPheno/pedmat.csv",stringsAsFactors=F,header=T)
allped$sire.pop[which(allped$sire.pop=="ml")]="m" 
allped$dam.pop[which(allped$dam.pop=="ml")]="m"
#this is so that they come out as being the 5th in the alphabet, because the pop that is malinalco in the provided coancestry matrix is row/column 5. 
#if it is re-run with population of malinalco 6th in genotype order and labeled as "ML", it would come out 6th in both coancestry and factor levels
allped$dam.pop<-as.numeric(as.factor(allped$dam.pop))
allped$sire.pop<- as.numeric(as.factor(allped$sire.pop))
colnames(allpheno)[1]<- "ID"
soils <- read.csv("~/GenoPheno/covmat.csv",stringsAsFactors=F,header=T)
#split data into three matrices by soil
traits.mt <-allpheno[which(soils$treatment==3),]
traits.mc <- allpheno[which(soils$treatment==5),]
traits.tc <- allpheno[which(soils$treatment==7),]
 
ped.mt <- allped[which(soils$treatment==3),] 
ped.mc <- allped[which(soils$treatment==5),]
ped.tc <- allped[which(soils$treatment==7),]

#create replicate traits, ped, matrices for each soil treatment, reset ID column to be original ID, not scaled ID.
traits.mt <- scale(traits.mt, center = TRUE, scale = TRUE)
traits.mc <- scale(traits.mc, center = TRUE, scale = TRUE)
traits.tc <- scale(traits.tc, center = TRUE, scale = TRUE)
traits.mt[,1] <-allpheno[which(soils$treatment==3),1]
traits.mc[,1] <- allpheno[which(soils$treatment==5),1]
traits.tc[,1] <- allpheno[which(soils$treatment==7),1]

#highland only version, subset, then center and scale, and reset ID column so it's not scaled
traits.mt.h <- traits.mt[which(ped.mt$dam.pop!=5 & ped.mt$dam.pop!=10),]
traits.mc.h <- traits.mc[which(ped.mc$dam.pop!=5 & ped.mc$dam.pop!=10),]
traits.tc.h <- traits.tc[which(ped.tc$dam.pop!=5 & ped.tc$dam.pop!=10),]
traits.mt.h <- scale(traits.mt.h, center = TRUE, scale = TRUE)
traits.mc.h <- scale(traits.mc.h, center = TRUE, scale = TRUE)
traits.tc.h <- scale(traits.tc.h, center = TRUE, scale = TRUE)
traits.mt.h[,1] <- traits.mt[which(ped.mt$dam.pop!=5 & ped.mt$dam.pop!=10),1]
traits.mc.h[,1] <- traits.mc[which(ped.mc$dam.pop!=5 & ped.mc$dam.pop!=10),1]
traits.tc.h[,1] <- traits.tc[which(ped.tc$dam.pop!=5 & ped.tc$dam.pop!=10),1]

#subset pedigree files for highland only analysis
ped.mt.h <- ped.mt[which(ped.mt$dam.pop!=5 & ped.mt$dam.pop!=10),] 
ped.mc.h <- ped.mc[which(ped.mc$dam.pop!=5 & ped.mc$dam.pop!=10),]
ped.tc.h <- ped.tc[which(ped.tc$dam.pop!=5 & ped.tc$dam.pop!=10),]
ped.mc.h$dam.pop[which(ped.mc.h$dam.pop>4)]<-ped.mc.h$dam.pop[which(ped.mc.h$dam.pop>4)]-1
ped.mc.h$sire.pop[which(ped.mc.h$sire.pop>4)]<-ped.mc.h$sire.pop[which(ped.mc.h$sire.pop>4)]-1
ped.mt.h$dam.pop[which(ped.mt.h$dam.pop>4)]<-ped.mt.h$dam.pop[which(ped.mt.h$dam.pop>4)]-1
ped.mt.h$sire.pop[which(ped.mt.h$sire.pop>4)]<-ped.mt.h$sire.pop[which(ped.mt.h$sire.pop>4)]-1
ped.tc.h$dam.pop[which(ped.tc.h$dam.pop>4)]<-ped.tc.h$dam.pop[which(ped.tc.h$dam.pop>4)]-1
ped.tc.h$sire.pop[which(ped.tc.h$sire.pop>4)]<-ped.tc.h$sire.pop[which(ped.tc.h$sire.pop>4)]-1

#driftsel needs, pedigree, covariate, and trait data, MUST be matrices
#structure of commands looks like:
#samp <- MH(afm$theta, ped, covars, traits, 100, 40, 2, alt=T) 

#these each take a very long time to run. i.e. can be weeks. 
#Each MH command at the very least should be run as a separate job in a separate script (i.e. duplicate this script, comment out all but one, run script, repeat for each dataset)
samp.mc <- MH(all.coan$theta,ped.mc,covars=traits.mc[,1],traits.mc,440000,40000,2000,alt=T)
save(samp.mc,file="driftselobjects/scaledtraits.samp.mc_med_LG95.Rdata")
samp.mt <- MH(all.coan$theta,ped.mt,covars=traits.mt[,1],traits.mt,440000,40000,2000,alt=T)
save(samp.mt,file="driftselobjects/scaledtraits.samp.mt_med_LG95.Rdata")
samp.tc <- MH(all.coan$theta,ped.tc,covars=traits.tc[,1],traits.tc,440000,40000,2000,alt=T)
save(samp.tc,file="driftselobjects/scaledtraits.samp.tc_med_LG95.Rdata")
#highlandONLY
highlandcoan <- all.coan$theta[c(1:4,6:9),c(1:4,6:9),]#
samp.mc.h <- MH(highlandcoan,ped.mc.h,covars=traits.mc.h[,1],traits.mc.h,440000,40000,2000,alt=T)#thinning in paper was 1000, does not change answer
save(samp.mc.h,file="driftselobjects/scaledtraits.samp.mc.h_med_LG95_440k.Rdata")
samp.mt.h <- MH(highlandcoan,ped.mt.h,covars=traits.mt.h[,1],traits.mt.h,440000,40000,2000,alt=T)#thinning in paper was 1000, does not change answer
save(samp.mt.h,file="driftselobjects/scaledtraits.samp.mt.h_med_LG95_440k.Rdata")
samp.tc.h <- MH(highlandcoan,ped.tc.h,covars=traits.tc.h[,1],traits.tc.h,440000,40000,2000,alt=T)
save(samp.tc.h,file="driftselobjects/lessLLL.samp.tc.h_med_LG95.Rdata")


#neutral tests
load("driftselobjects/scaledtraits.samp.mt_med_LG95")
load("driftselobjects/scaledtraits.samp.mc_med_LG95.Rdata")
load("driftselobjects/scaledtraits.samp.tc_med_LG95.Rdata")
load("driftselobjects/scaledtraits.samp.mt.h_med_LG95_440k.Rdata")
load("driftselobjects/scaledtraits.samp.mc.h_med_LG95_440k.Rdata")
load("driftselobjects/scaledtraits.samp.tc.h_med_LG95_440k.Rdata")
nt.mc <- neut.test(samp.mc$pop.ef, samp.mc$G, samp.mc$theta, silent=T)
nt.mt <- neut.test(samp.mt$pop.ef, samp.mt$G, samp.mt$theta, silent=T)
nt.tc <- neut.test(samp.tc$pop.ef, samp.tc$G, samp.tc$theta, silent=T)
nt.mc.h <- neut.test(samp.mc.h$pop.ef, samp.mc.h$G, samp.mc.h$theta, silent=T)
nt.mt.h <- neut.test(samp.mt.h$pop.ef, samp.mt.h$G, samp.mt.h$theta, silent=T)
nt.tc.h <- neut.test(samp.tc.h$pop.ef, samp.tc.h$G, samp.tc.h$theta, silent=T)

results <- data.frame(label =c("MC","MT","TC"),test =  c(nt.mc,nt.mt,nt.tc))
sink("~/GenoPheno/scaled.results_mLG95.txt")
print(results)
sink()

results.h <- data.frame(label =c("MC","MT","TC"),test =  c(nt.mc.h,nt.mt.h,nt.tc.h))
sink("~/GenoPheno/scaled.h.results_HmLG95.txt")
print(results)
sink()

# H tests (environment tests)
envmatall <-read.csv("~/GenoPheno/AlphabeticalPopEnvDat.csv")
envmat <- envmatall[,c(1,4,48)]#just mean temp, annual precip, and soil water holding capacity


Hmt <-  H.test(samp.mt$pop.ef, samp.mt$G, samp.mt$theta, envmat, silent=TRUE)
Hmc <-  H.test(samp.mc$pop.ef, samp.mc$G, samp.mc$theta, envmat, silent=TRUE)
Htc <-  H.test(samp.tc$pop.ef, samp.tc$G, samp.tc$theta, envmat, silent=TRUE)

Hmt.h <-  H.test(samp.mt.h$pop.ef, samp.mt.h$G, samp.mt.h$theta, envmat[-c(5,10),], silent=TRUE)
Hmc.h <-  H.test(samp.mc.h$pop.ef, samp.mc.h$G, samp.mc.h$theta, envmat[-c(5,10),], silent=TRUE)
Htc.h <-  H.test(samp.tc.h$pop.ef, samp.tc.h$G, samp.tc.h$theta, envmat[-c(5,10),], silent=TRUE)

resultH <- data.frame(names=c("Hmt","Hmc","Htc"),Htest.result = c(Hmt,Hmc,Htc))
resultH.h <- data.frame(names=c("Hmt.h","Hmc.h","Htc.h"),Htest.result = c(Hmt.h,Hmc.h,Htc.h))

sink("~/GenoPheno/resultH_med_LG95.txt")
print(resultH)
sink()

sink("~/GenoPheno/resultH_med_LG95.h.txt")
print(resultH.h)
sink()

