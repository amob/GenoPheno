##Fit models with random effects of mom, population, and mom and population in soil treatment
##Testing also whether adding fixed effects of environmental variation explains population effects 

library(MCMCglmm)

range01=function(x){
newnums <- sapply(1:length(x), function(z) ifelse(!is.na(x[z]), (x[z]-min(x,na.rm=T))/(max(x,na.rm=T)-min(x,na.rm=T)), NA ))
return(newnums)
}

fitmodelset <- function(traitcols,mainvarformula,dataframe) { #names of trait variables in data frame, fixed effects formula in text e.g. "~ 1", "~elevation", and data
	rbase<- list()
	rpxs <- list()
	rfxs <- list()
	rfxspxs <- list()
	randbase <- as.formula("~  dam.pop + dam + soiltrt")
	randpxs <- as.formula("~  dam.pop + dam + soiltrt + dam.pop:soiltrt") 
	randfxs <- as.formula("~  dam.pop + dam + soiltrt + dam:soiltrt") 
	randfxspxs <- as.formula("~  dam.pop + dam + soiltrt + dam:soiltrt + dam.pop:soiltrt") 
	for(i in 1:length(traitcols)){
		mainform <- as.formula(paste(traitcols[i], mainvarformula,sep=""))
		rbase[[i]] <- lapply(1:10, function(z)  
			MCMCglmm(mainform, random=randbase, family="gaussian", data=alldat,verbose=FALSE,pr=TRUE,nitt=100000, thin=100, burnin=20000))
		rpxs[[i]] <- lapply(1:10, function(z)  
			MCMCglmm(mainform, random=randpxs, family="gaussian", data=alldat,verbose=FALSE,pr=TRUE,nitt=100000, thin=100, burnin=20000))
		rfxs[[i]] <- lapply(1:10, function(z)  
			MCMCglmm(mainform, random=randfxs, family="gaussian", data=alldat,verbose=FALSE,pr=TRUE,nitt=100000, thin=100, burnin=20000))
		rfxspxs[[i]] <- lapply(1:10, function(z)  
			MCMCglmm(mainform, random=randfxspxs, family="gaussian",data=alldat,verbose=FALSE,pr=TRUE,nitt=100000, thin=100, burnin=20000))
	}
	clusterlist <- list(rbase,rpxs,rfxs,rfxspxs)
	return(clusterlist)
}

mod.summs <- function(filepath,ntraits,nmodtype,nreps) { #lazy issue in that filepath is to saved rdata file with "clusterlist" object, structured exactly as the environmental ones
	load(filepath)
	diccl <- c()
	pval <- c()
	slope <- c()
	counter <- 1
	for(k in 1:ntraits){#traits inside lists of model type
		for(j in 1:nmodtype){#model type is the first level of list
			for(i in 1:nreps){#repeats innermost level
				diccl[counter] <- clusterlist[[j]][[k]][[i]]$DIC
				pval[counter] <- summary(clusterlist[[j]][[k]][[i]])$solutions[2,5]
				slope[counter] <- summary(clusterlist[[j]][[k]][[i]])$solutions[2,1]
				counter <- counter +1
				}
			}
		}
	arraydiccl<-array(diccl,dim=c(10,4,9))#fills first array first(last dimension), first row first(first dimension)
	arraypval <-array(pval,dim=c(10,4,9))
	arrayslp <-array(slope,dim=c(10,4,9))
	rm(clusterlist)
	return(list(arraydiccl,arraypval,arrayslp))
}

####READ IN DATA
soils <- read.csv("~/GenoPheno/covmat.csv",stringsAsFactors=F,header=T)
allpheno <- read.csv("~/GenoPheno/phenomat.csv",stringsAsFactors=F,header=T)
allped <- read.csv("~/GenoPheno/pedmat.csv",stringsAsFactors=F,header=T)
Apopenvdat <-read.csv("AlphabeticalPopEnvDat.csv",stringsAsFactors=F,header=T)

greens <- (allpheno$green)/(allpheno$red + allpheno$blue + allpheno$green)
#swaps colors to ratios with total brightness
allpheno <- allpheno[,-c(3,10:12,14,16)]#cuts problematic traits
allpheno$green <- greens
colnames(allpheno)[1]<- "ID"
cpheno <- scale(allpheno[,-1])
#soil treatment varible
soiltrt <- c()
soiltrt[which(soils$treatment==3)] <- "mt"
soiltrt[which(soils$treatment==5)] <- "mc"
soiltrt[which(soils$treatment==7)] <- "tc"
#soil trt and ped numbers reflect order in which populations were SAMPLED, not alphabet.

allped$sire.pop[which(allped$sire.pop=="ml")]="m"
allped$dam.pop[which(allped$dam.pop=="ml")]="m"#"m"is 5th in the factorized, malinalco in the coancestry matrix is row/column 5. "ML" would come out 6th.
allped$dam.pop<-as.numeric(as.factor(allped$dam.pop))
allped$sire.pop<- as.numeric(as.factor(allped$sire.pop))
#now ped pop numbers reflect alphabet, although dams still reflect field sample order, and sire #s are unrelated and unique for all

envdat <- Apopenvdat[allped$dam.pop,]
cenvdat <- scale(envdat)
#since dam.pop now reflects alphabetical order, can pull out env info this way

#### FIGURE
## commented out because done, no need to rewrite file
##ploting pop-momxe effects.
#plot phenotypes across soil treatments
# color by pop source temp. plot all individuals connected by lines?
traits.mt <-allpheno[which(soils$treatment==3),]
traits.mc <- allpheno[which(soils$treatment==5),]
traits.tc <- allpheno[which(soils$treatment==7),]
ped.mt <- allped[which(soils$treatment==3),] 
ped.mc <- allped[which(soils$treatment==5),]
ped.tc <- allped[which(soils$treatment==7),]
damtrait.mc <- sapply(2:10,function(z) tapply(traits.mc[,z],ped.mc$dam,mean,na.rm=T))
damtrait.mt <- sapply(2:10,function(z) tapply(traits.mt[,z],ped.mt$dam,mean,na.rm=T))
damtrait.tc <- sapply(2:10,function(z) tapply(traits.tc[,z],ped.tc$dam,mean,na.rm=T))
poptrait.mc <- sapply(2:10,function(z) tapply(traits.mc[,z],ped.mc$dam.pop,mean,na.rm=T))
poptrait.mt <- sapply(2:10,function(z) tapply(traits.mt[,z],ped.mt$dam.pop,mean,na.rm=T))
poptrait.tc <- sapply(2:10,function(z) tapply(traits.tc[,z],ped.tc$dam.pop,mean,na.rm=T))
damtraits.soils <- array(,dim=c(100,3,9))#100 moms, 3 soils, 9 traits ##REORDERED SOILS
for(i in 1:9){
damtraits.soils[,2,i] <- damtrait.mc[,i]
damtraits.soils[,1,i] <- damtrait.mt[,i]
damtraits.soils[,3,i] <- damtrait.tc[,i]
}
poptraits.soils <- array(,dim=c(10,3,9))#10 pops, 3 soils, 9 traits ## REORDERED SOILS
for(i in 1:9){
poptraits.soils[,2,i] <- poptrait.mc[,i]
poptraits.soils[,1,i] <- poptrait.mt[,i]
poptraits.soils[,3,i] <- poptrait.tc[,i]
}
traitcols <- c("Days to flowering","Days to germination","Tassel length cm","Shoot biomass g","Root biomass g","Height cm","Stem width mm","Leaf width mm","Stem greenness %")#nice names
rbb<-rgb(range01(Apopenvdat$TAnn),0,1-range01(Apopenvdat$TAnn))
rbsoft2<- rgb(range01(Apopenvdat$TAnn),0,1-range01(Apopenvdat$TAnn),alpha=.1)
pdf("~/GenoPheno/popdamtraits_soilsOct2018.pdf",height=11,width=11)
layout(matrix(1:12,ncol=3,byrow=T),heights=c(5,5,5,2))
par(mar=c(4,6,1,2))
for(i in 1:6){
matplot(t(damtraits.soils[,,i]),type="l",col=rep(rbsoft2)[damtoalphpop],lty=1,ylab=traitcols[i],xlab="",xaxt="n",cex.axis=2.25,cex.lab=2.75)
matplot(t(poptraits.soils[,,i]),type="l",col=rbb,lty=1,add=T)
}
for(i in 7:9){
matplot(t(damtraits.soils[,,i]),type="l",col=rep(rbsoft2)[damtoalphpop],lty=1,ylab=traitcols[i],xlab="",xaxt="n",cex.axis=2.25,cex.lab=2.75)
matplot(t(poptraits.soils[,,i]),type="l",col=rbb,lty=1,add=T)
axis(side=1,at=c(1,2,3),labels=c("Biota15.0","Biota14.3","Biota13.0"),cex.axis=2.75,las=2)
}
dev.off()

 

# #ANOVA analysis, finds GxE at family, not pop level explains sig variation for many traits.
# #allofit <- cbind(ctrait,cped,soiltrt)
sse.dam <- list()
for(i in 1:9){
sse.dam[[i]] <- summary(aov(cpheno[,i]~as.character(allped$dam) +soiltrt + as.character(allped$dam):soiltrt,,na.action=na.omit))#only allows dam or pop at once
}#none of the pop ints and lots (about half of traits) of dam ints (abv) are sig -- lots of GxE
sse.pop <- list()
for(i in 1:9){
sse.pop[[i]] <- summary(aov(cpheno[,i]~as.character(allped$dam.pop) +soiltrt + as.character(allped$dam.pop):soiltrt,na.action=na.omit))
}#
# 


####MODELS

alldat <- cbind(cpheno,allped,soiltrt,cenvdat)
alldat$dam <- as.character(alldat$dam)
alldat$dam.pop <- as.character(alldat$dam.pop)
traitcols <- c("flowerday","germday","espigalength","shoot","root","finalheight","finalstemw","width3leaf","green")#functional names, matching phenotype data frame

##INTERCEPT
clusterlist <- fitmodelset(traitcols,"~1",alldat)
save(clusterlist,file="~/GenoPheno/traitgxemodelsINTOct2018-cL.Rdata")
# ELEVATION
clusterlist <- fitmodelset(traitcols,"~ elevation",alldat)
save(clusterlist,file="~/GenoPheno/traitgxemodelsELEVOct2018-cL.Rdata")
#MAT
clusterlist <- fitmodelset(traitcols,"~ TAnn",alldat)
save(clusterlist,file="~/GenoPheno/traitgxemodelsTAnnOct2018-cL.Rdata")
# Pann
clusterlist <- fitmodelset(traitcols,"~ Pann",alldat)
save(clusterlist,file="~/GenoPheno/traitgxemodelsPannOct2018-cL.Rdata")
#SWC 
clusterlist <- fitmodelset(traitcols,"~ CapacidadDeCampoPercent",alldat)
save(clusterlist,file="~/GenoPheno/traitgxemodelsSWCOct2018-cL.Rdata")

load("~/GenoPheno/traitgxemodelsINTOct2018-cL.Rdata")
diccl <- c()
counter <- 1
for(k in 1:9){#traits inside lists of model type
for(j in 1:4){#model type is the first level of list
for(i in 1:10){#repeats innermost level
	diccl[counter] <- clusterlist[[j]][[k]][[i]]$DIC
	counter <- counter +1
	}
}
}
rm(clusterlist)
arraydicclI<-array(diccl,dim=c(10,4,9))#fills first array first(last dimension), first row first(first dimension)

MATmods <- mod.summs("~/GenoPheno/traitgxemodelsTAnnOct2018-cL.Rdata",9,4,10)
TAPmods <- mod.summs("~/GenoPheno/traitgxemodelsPannOct2018-cL.Rdata",9,4,10)
SWCmods <- mod.summs("~/GenoPheno/traitgxemodelsSWCOct2018-cL.Rdata",9,4,10)
ELVmods <- mod.summs("~/GenoPheno/traitgxemodelsELEVOct2018-cL.Rdata",9,4,10)


# ##SUMMARY
diclist <- list( apply(arraydicclI,c(2,3),mean), apply(MATmods[[1]],c(2,3),mean), apply(TAPmods[[1]],c(2,3),mean),
		apply(SWCmods[[1]],c(2,3),mean), apply(ELVmods[[1]],c(2,3),mean) )

pvallist <- list( apply(MATmods[[2]],c(2,3),mean), apply(TAPmods[[2]],c(2,3),mean),
	  apply(SWCmods[[2]],c(2,3),mean), apply(ELVmods[[2]],c(2,3),mean) )

slopelist <- list( apply(MATmods[[3]],c(2,3),mean), apply(TAPmods[[3]],c(2,3),mean),
	  apply(SWCmods[[3]],c(2,3),mean), apply(ELVmods[[3]],c(2,3),mean) )


resultlist <- list(diclist,pvallist,slopelist)
names(resultlist) <- c("dic","pvals","slopes")
sink("~/GenoPheno/04_02_gxeresultsOct2018-cL.txt")
print(resultlist)
sink()
print(resultlist)
