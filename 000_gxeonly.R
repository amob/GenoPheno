##Fit models with random effects of mom, population, and mom and population in soil treatment
##Testing also whether adding fixed effects of environmental variation explains population effects 

library(MCMCglmm)

range01=function(x){
newnums <- sapply(1:length(x), function(z) ifelse(!is.na(x[z]), (x[z]-min(x,na.rm=T))/(max(x,na.rm=T)-min(x,na.rm=T)), NA ))
return(newnums)
}

soils <- read.csv("~/GenoPheno/covmat.csv",stringsAsFactors=F,header=T)
allpheno <- read.csv("~/GenoPheno/phenomat.csv",stringsAsFactors=F,header=T)
greens <- (allpheno$green)/(allpheno$red + allpheno$blue + allpheno$green) #swaps color to ratio with total brightness
allpheno <- allpheno[,-c(3,10:12,14,16)]#cuts traits
allpheno$green <- greens
colnames(allpheno)[1]<- "ID"


allped <- read.csv("~/GenoPheno/pedmat.csv",stringsAsFactors=F,header=T)
allped$sire.pop[which(allped$sire.pop=="ml")]="m"
allped$dam.pop[which(allped$dam.pop=="ml")]="m"#"m"is 5th in the factorized, malinalco in the coancestry matrix is row/column 5. "ML" would come out 6th.
allped$dam.pop<-as.numeric(as.factor(allped$dam.pop))
allped$sire.pop<- as.numeric(as.factor(allped$sire.pop))

Apopenvdat <-read.csv("AlphabeticalPopEnvDat.csv",stringsAsFactors=F,header=T)
envdat <- Apopenvdat[allped$dam.pop,]

#this just reorganizes data in a way that makes more sense
traits.mt <-allpheno[which(soils$treatment==3),]
traits.mc <- allpheno[which(soils$treatment==5),]
traits.tc <- allpheno[which(soils$treatment==7),]
ped.mt <- allped[which(soils$treatment==3),] 
env.mt <- envdat[which(soils$treatment==3),] 
ped.mc <- allped[which(soils$treatment==5),]
env.mc <- envdat[which(soils$treatment==5),] 
ped.tc <- allped[which(soils$treatment==7),]
env.tc <- envdat[which(soils$treatment==7),] 
#listpops<- ("CL","CU","DA","FP","M","MC","MT","TC","TX","TZ")#"M" = "ML"

ctrait <- rbind(scale(traits.mc[,2:10]),scale(traits.mt[,2:10]),scale(traits.tc[,2:10]))
cped <- rbind(ped.mc,ped.mt,ped.tc)
soiltrt <- c(rep("mc",times=nrow(ped.mc)),rep("mt",times=nrow(ped.mt)),rep("tc",times=nrow(ped.tc)))
ordenv <- rbind(env.mc,env.mt,env.tc)

alldat <- cbind(ctrait,cped,soiltrt,ordenv)
alldat$dam <- as.character(alldat$dam)
alldat$dam.pop <- as.character(alldat$dam.pop)

########
##FIGURE
damtrait.mc <- sapply(2:10,function(z) tapply(traits.mc[,z],ped.mc$dam,mean,na.rm=T))
damtrait.mt <- sapply(2:10,function(z) tapply(traits.mt[,z],ped.mt$dam,mean,na.rm=T))
damtrait.tc <- sapply(2:10,function(z) tapply(traits.tc[,z],ped.tc$dam,mean,na.rm=T))

poptrait.mc <- sapply(2:10,function(z) tapply(traits.mc[,z],ped.mc$dam.pop,mean,na.rm=T))
poptrait.mt <- sapply(2:10,function(z) tapply(traits.mt[,z],ped.mt$dam.pop,mean,na.rm=T))
poptrait.tc <- sapply(2:10,function(z) tapply(traits.tc[,z],ped.tc$dam.pop,mean,na.rm=T))

damtraits.soils <- array(,dim=c(100,3,9))#100 moms, 3 soils, 9 traits ##REORDERED SOILS so mt (warmest source) first
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

traitcols <- c("Days to flowering","Days to germination","Tassel length cm","Shoot biomass g","Root biomass g","Height cm","Stem width mm","Leaf width mm","Stem greenness %")
rbb<-rgb(range01(Apopenvdat$TAnn),0,1-range01(Apopenvdat$TAnn))

rbsoft2<- rgb(range01(Apopenvdat$TAnn),0,1-range01(Apopenvdat$TAnn),alpha=.1)
pdf("~/GenoPheno/popdamtraits_soils.pdf",height=11,width=11)
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

##MODELS FOR TABLE 
#the following are the repeated similar sets of models. A function could be written to shorten this code.

#RAND ONLY
traitcols <- c("flowerday","germday","espigalength","shoot","root","finalheight","finalstemw","width3leaf","green")
baselist <- list()
poplist <- list()
damlist <- list()
dampoplist <- list()
for(i in 1:length(traitcols)){
	mainform <-as.formula(paste(traitcols[i],"~ soiltrt",sep=""))
	randbase <- as.formula("~  dam.pop:dam")
	randpop <- as.formula("~   dam.pop:dam + dam.pop:soiltrt")
	randdam <- as.formula("~   dam.pop:dam:soiltrt")
	randdampop <- as.formula("~dam.pop:dam:soiltrt + dam.pop:soiltrt")
	baselist[[i]] <- lapply(1:10, function(z) MCMCglmm(mainform, random= randbase, family="gaussian",
		data=alldat,verbose=FALSE,pr=TRUE,nitt=100000, thin=50, burnin=25000))
	poplist[[i]] <- lapply(1:10, function(z) MCMCglmm(mainform, random= randpop, family="gaussian",
		data=alldat,verbose=FALSE,pr=TRUE,nitt=100000, thin=50, burnin=25000))
	damlist[[i]] <- lapply(1:10, function(z) MCMCglmm(mainform, random= randdam, family="gaussian",
		data=alldat,verbose=FALSE,pr=TRUE,nitt=100000, thin=50, burnin=25000))
	dampoplist[[i]] <- lapply(1:10, function(z) MCMCglmm(mainform, random= randdampop, family="gaussian",
		data=alldat,verbose=FALSE,pr=TRUE,nitt=100000, thin=50, burnin=25000))
}
clusterlist <- list(baselist,damlist,poplist,dampoplist)
save(clusterlist,file="~/GenoPheno/traitgxemodels.Rdata")
load("~/GenoPheno/traitgxemodels.Rdata")
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
arraydiccl<-array(diccl,dim=c(10,4,9))#fills first array first(last dimension), first row first(first dimension)
rm(diccl)
rm(clusterlist)

##WITH FIXED POP ELEVATION
traitcols <- c("flowerday","germday","espigalength","shoot","root","finalheight","finalstemw","width3leaf","green")
baselist <- list()
poplist <- list()
damlist <- list()
dampoplist <- list()
for(i in 1:length(traitcols)){
	mainform <-as.formula(paste(traitcols[i],"~ soiltrt + elevation",sep=""))
	randbase <- as.formula("~  dam.pop:dam")
	randpop <- as.formula("~   dam.pop:dam + dam.pop:soiltrt")
	randdam <- as.formula("~   dam.pop:dam:soiltrt")
	randdampop <- as.formula("~dam.pop:dam:soiltrt + dam.pop:soiltrt")
	baselist[[i]] <- lapply(1:10, function(z) MCMCglmm(mainform, random= randbase, family="gaussian",
		data=alldat,verbose=FALSE,pr=TRUE,nitt=100000, thin=50, burnin=25000))
	poplist[[i]] <- lapply(1:10, function(z) MCMCglmm(mainform, random= randpop, family="gaussian",
		data=alldat,verbose=FALSE,pr=TRUE,nitt=100000, thin=50, burnin=25000))
	damlist[[i]] <- lapply(1:10, function(z) MCMCglmm(mainform, random= randdam, family="gaussian",
		data=alldat,verbose=FALSE,pr=TRUE,nitt=100000, thin=50, burnin=25000))
	dampoplist[[i]] <- lapply(1:10, function(z) MCMCglmm(mainform, random= randdampop, family="gaussian",
		data=alldat,verbose=FALSE,pr=TRUE,nitt=100000, thin=50, burnin=25000))
}
clusterlist <- list(baselist,damlist,poplist,dampoplist)
save(clusterlist,file="~/GenoPheno/traitgxemodelsELEV.Rdata")
load("~/GenoPheno/traitgxemodelsELEV.Rdata")
diccl <- c()
pvalEnv <- c()
slopeEnv <- c()
counter <- 1
for(k in 1:9){#traits inside lists of model type
for(j in 1:4){#model type is the first level of list
for(i in 1:10){#repeats innermost level
	diccl[counter] <- clusterlist[[j]][[k]][[i]]$DIC
	pvalEnv[counter] <- summary(clusterlist[[j]][[k]][[i]])$solutions[4,5]
	slopeEnv[counter] <- summary(clusterlist[[j]][[k]][[i]])$solutions[4,1]
	counter <- counter +1
	}
}
}
arraydicclE <-array(diccl,dim=c(10,4,9))#fills first array first(last dimension), first row first(first dimension)
arraypvalE <-array(pvalEnv,dim=c(10,4,9))
arrayslpE <-array(slopeEnv,dim=c(10,4,9))
rm(diccl); rm(pvalEnv); rm(slopeEnv)
rm(clusterlist)

#WITH FIXED POP MAT
traitcols <- c("flowerday","germday","espigalength","shoot","root","finalheight","finalstemw","width3leaf","green")
baselist <- list()
poplist <- list()
damlist <- list()
dampoplist <- list()
for(i in 1:length(traitcols)){
	mainform <-as.formula(paste(traitcols[i],"~ soiltrt + TAnn",sep=""))
	randbase <- as.formula("~  dam.pop:dam")
	randpop <- as.formula("~   dam.pop:dam + dam.pop:soiltrt")
	randdam <- as.formula("~   dam.pop:dam:soiltrt")
	randdampop <- as.formula("~dam.pop:dam:soiltrt + dam.pop:soiltrt")
	baselist[[i]] <- lapply(1:10, function(z) MCMCglmm(mainform, random= randbase, family="gaussian",
		data=alldat,verbose=FALSE,pr=TRUE,nitt=100000, thin=50, burnin=25000))
	poplist[[i]] <- lapply(1:10, function(z) MCMCglmm(mainform, random= randpop, family="gaussian",
		data=alldat,verbose=FALSE,pr=TRUE,nitt=100000, thin=50, burnin=25000))
	damlist[[i]] <- lapply(1:10, function(z) MCMCglmm(mainform, random= randdam, family="gaussian",
		data=alldat,verbose=FALSE,pr=TRUE,nitt=100000, thin=50, burnin=25000))
	dampoplist[[i]] <- lapply(1:10, function(z) MCMCglmm(mainform, random= randdampop, family="gaussian",
		data=alldat,verbose=FALSE,pr=TRUE,nitt=100000, thin=50, burnin=25000))
}
clusterlist <- list(baselist,damlist,poplist,dampoplist)
save(clusterlist,file="~/GenoPheno/traitgxemodelsTAnn.Rdata")
load("~/GenoPheno/traitgxemodelsTAnn.Rdata")
diccl <- c()
pvalEnv <- c()
slopeEnv <- c()
counter <- 1
for(k in 1:9){#traits inside lists of model type
for(j in 1:4){#model type is the first level of list
for(i in 1:10){#repeats innermost level
	diccl[counter] <- clusterlist[[j]][[k]][[i]]$DIC
	pvalEnv[counter] <- summary(clusterlist[[j]][[k]][[i]])$solutions[4,5]
	slopeEnv[counter] <- summary(clusterlist[[j]][[k]][[i]])$solutions[4,1]
	counter <- counter +1
	}
}
}
arraydicclT <-array(diccl,dim=c(10,4,9))#fills first array first(last dimension), first row first(first dimension)
arraypvalT <-array(pvalEnv,dim=c(10,4,9))
arrayslpT <-array(slopeEnv,dim=c(10,4,9))
rm(diccl); rm(pvalEnv); rm(slopeEnv)
rm(clusterlist)

##WITH FIXED POP Pann
traitcols <- c("flowerday","germday","espigalength","shoot","root","finalheight","finalstemw","width3leaf","green")
baselist <- list()
poplist <- list()
damlist <- list()
dampoplist <- list()
for(i in 1:length(traitcols)){
	mainform <-as.formula(paste(traitcols[i],"~ soiltrt + Pann",sep=""))
	randbase <- as.formula("~  dam.pop:dam")
	randpop <- as.formula("~   dam.pop:dam + dam.pop:soiltrt")
	randdam <- as.formula("~   dam.pop:dam:soiltrt")
	randdampop <- as.formula("~dam.pop:dam:soiltrt + dam.pop:soiltrt")
	baselist[[i]] <- lapply(1:10, function(z) MCMCglmm(mainform, random= randbase, family="gaussian",
		data=alldat,verbose=FALSE,pr=TRUE,nitt=100000, thin=50, burnin=25000))
	poplist[[i]] <- lapply(1:10, function(z) MCMCglmm(mainform, random= randpop, family="gaussian",
		data=alldat,verbose=FALSE,pr=TRUE,nitt=100000, thin=50, burnin=25000))
	damlist[[i]] <- lapply(1:10, function(z) MCMCglmm(mainform, random= randdam, family="gaussian",
		data=alldat,verbose=FALSE,pr=TRUE,nitt=100000, thin=50, burnin=25000))
	dampoplist[[i]] <- lapply(1:10, function(z) MCMCglmm(mainform, random= randdampop, family="gaussian",
		data=alldat,verbose=FALSE,pr=TRUE,nitt=100000, thin=50, burnin=25000))
}
clusterlist <- list(baselist,damlist,poplist,dampoplist)
save(clusterlist,file="~/GenoPheno/traitgxemodelsPann.Rdata")
load("~/GenoPheno/traitgxemodelsPann.Rdata")
diccl <- c()
pvalEnv <- c()
slopeEnv <- c()
counter <- 1
for(k in 1:9){#traits inside lists of model type
for(j in 1:4){#model type is the first level of list
for(i in 1:10){#repeats innermost level
	diccl[counter] <- clusterlist[[j]][[k]][[i]]$DIC
	pvalEnv[counter] <- summary(clusterlist[[j]][[k]][[i]])$solutions[4,5]
	slopeEnv[counter] <- summary(clusterlist[[j]][[k]][[i]])$solutions[4,1]
	counter <- counter +1
	}
}
}
arraydicclP <-array(diccl,dim=c(10,4,9))#fills first array first(last dimension), first row first(first dimension)
arraypvalP <-array(pvalEnv,dim=c(10,4,9))
arrayslpP <-array(slopeEnv,dim=c(10,4,9))
rm(diccl); rm(pvalEnv); rm(slopeEnv)
rm(clusterlist)

##WITH FIXED POP SWC
traitcols <- c("flowerday","germday","espigalength","shoot","root","finalheight","finalstemw","width3leaf","green")
baselist <- list()
poplist <- list()
damlist <- list()
dampoplist <- list()
for(i in 1:length(traitcols)){
	mainform <-as.formula(paste(traitcols[i],"~ soiltrt + CapacidadDeCampoPercent",sep=""))
	randbase <- as.formula("~  dam.pop:dam")
	randpop <- as.formula("~   dam.pop:dam + dam.pop:soiltrt")
	randdam <- as.formula("~   dam.pop:dam:soiltrt")
	randdampop <- as.formula("~dam.pop:dam:soiltrt + dam.pop:soiltrt")
	baselist[[i]] <- lapply(1:10, function(z) MCMCglmm(mainform, random= randbase, family="gaussian",
		data=alldat,verbose=FALSE,pr=TRUE,nitt=100000, thin=50, burnin=25000))
	poplist[[i]] <- lapply(1:10, function(z) MCMCglmm(mainform, random= randpop, family="gaussian",
		data=alldat,verbose=FALSE,pr=TRUE,nitt=100000, thin=50, burnin=25000))
	damlist[[i]] <- lapply(1:10, function(z) MCMCglmm(mainform, random= randdam, family="gaussian",
		data=alldat,verbose=FALSE,pr=TRUE,nitt=100000, thin=50, burnin=25000))
	dampoplist[[i]] <- lapply(1:10, function(z) MCMCglmm(mainform, random= randdampop, family="gaussian",
		data=alldat,verbose=FALSE,pr=TRUE,nitt=100000, thin=50, burnin=25000))
}
clusterlist <- list(baselist,damlist,poplist,dampoplist)
save(clusterlist,file="~/GenoPheno/traitgxemodelsSWC.Rdata")
load("~/GenoPheno/traitgxemodelsSWC.Rdata")
diccl <- c()
pvalEnv <- c()
slopeEnv <- c()
counter <- 1
for(k in 1:9){#traits inside lists of model type
for(j in 1:4){#model type is the first level of list
for(i in 1:10){#repeats innermost level
	diccl[counter] <- clusterlist[[j]][[k]][[i]]$DIC
	pvalEnv[counter] <- summary(clusterlist[[j]][[k]][[i]])$solutions[4,5]
	slopeEnv[counter] <- summary(clusterlist[[j]][[k]][[i]])$solutions[4,1]
	counter <- counter +1
	}
}
}
arraydicclSWC <-array(diccl,dim=c(10,4,9))#fills first array first(last dimension), first row first(first dimension)
arraypvalSWC <-array(pvalEnv,dim=c(10,4,9))
arrayslpSWC <-array(slopeEnv,dim=c(10,4,9))
rm(diccl); rm(pvalEnv); rm(slopeEnv)
rm(clusterlist)
 
# ##SUMMARY
diclist <- list( apply(arraydiccl,c(2,3),mean), apply(arraydicclT,c(2,3),mean), apply(arraydicclP,c(2,3),mean),
		apply(arraydicclSWC,c(2,3),mean), apply(arraydicclE,c(2,3),mean) )

pvallist <- list( apply(arraypvalT,c(2,3),mean), apply(arraypvalP,c(2,3),mean),
	  apply(arraypvalSWC,c(2,3),mean), apply(arraypvalE,c(2,3),mean) )

slopelist <- list( apply(arrayslpT,c(2,3),mean), apply(arrayslpP,c(2,3),mean),
	  apply(arrayslpSWC,c(2,3),mean), apply(arrayslpE,c(2,3),mean) )
resultlist <- list(diclist,pvallist,slopelist)
names(resultlist) <- c("dic","pvals","slopes")
sink("~/GenoPheno/04_02_gxeresults.txt")
print(resultlist)
sink()
print(resultlist)
