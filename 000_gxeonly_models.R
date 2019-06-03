##Fit models with random effects of mom, population, and mom and population in soil treatment
##Testing also whether adding fixed effects of environmental variation explains population effects 
library(MCMCglmm)

range01=function(x){
newnums <- sapply(1:length(x), function(z) ifelse(!is.na(x[z]), (x[z]-min(x,na.rm=T))/(max(x,na.rm=T)-min(x,na.rm=T)), NA ))
return(newnums)
}

fitmodelset <- function(traitcols,mainvarformula,dataframe) { #names of trait variables in data frame, fixed effects formula in text e.g. "~ 1", "~elevation", and data
	r2fp <- list()
	rbase <- list()
	rpxs <- list()
	rfxs <- list()
	rfxspxs <- list()
	rand2fp <- as.formula("~  dam + dam.pop")
	randbase <- as.formula("~  dam.pop + dam + soiltrt")
	randpxs <- as.formula("~  dam.pop + dam + soiltrt + dam.pop:soiltrt") 
	randfxs <- as.formula("~  dam.pop + dam + soiltrt + dam:soiltrt") 
	randfxspxs <- as.formula("~  dam.pop + dam + soiltrt + dam:soiltrt + dam.pop:soiltrt") 
	for(i in 1:length(traitcols)){
		mainform <- as.formula(paste(traitcols[i], mainvarformula,sep=""))
		r2fp[[i]] <- lapply(1:10, function(z)  
			MCMCglmm(mainform, random=rand2fp, family="gaussian", data=dataframe,verbose=FALSE,pr=TRUE,nitt=100000, thin=100, burnin=20000))
		rbase[[i]] <- lapply(1:10, function(z)  
			MCMCglmm(mainform, random=randbase, family="gaussian", data=dataframe,verbose=FALSE,pr=TRUE,nitt=100000, thin=100, burnin=20000))
		rpxs[[i]] <- lapply(1:10, function(z)  
			MCMCglmm(mainform, random=randpxs, family="gaussian", data=dataframe,verbose=FALSE,pr=TRUE,nitt=100000, thin=100, burnin=20000))
		rfxs[[i]] <- lapply(1:10, function(z)  
			MCMCglmm(mainform, random=randfxs, family="gaussian", data=dataframe,verbose=FALSE,pr=TRUE,nitt=100000, thin=100, burnin=20000))
		rfxspxs[[i]] <- lapply(1:10, function(z)  
			MCMCglmm(mainform, random=randfxspxs, family="gaussian",data=dataframe,verbose=FALSE,pr=TRUE,nitt=100000, thin=100, burnin=20000))
	}
	clusterlist <- list(r2fp,rbase,rpxs,rfxs,rfxspxs)
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
	arraydiccl<-array(diccl,dim=c(nreps,nmodtype,ntraits))#fills first array first(last dimension), first row first(first dimension)
	arraypval <-array(pval,dim=c(nreps,nmodtype,ntraits))
	arrayslp <-array(slope,dim=c(nreps,nmodtype,ntraits))
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
allpheno <- allpheno[,-c(3,10:12,14,16)]#cuts some traits (some colinear, some non-normal, some invariant), leaves 9 traits
allpheno$green <- greens
allpheno$germday <- allpheno$germday - 16 # is in days since july 1st, and later want it in Days to Germination, so subtract planting date

colnames(allpheno)[1]<- "ID"
cpheno <- scale(allpheno[,-1]) #center and scale
#soil treatment varible
soiltrt <- c()
soiltrt[which(soils$treatment==3)] <- "mt"
soiltrt[which(soils$treatment==5)] <- "mc"
soiltrt[which(soils$treatment==7)] <- "tc"
#soil trt and ped numbers reflect order in which populations were SAMPLED, not alphabet.

#redefining levels in pedigree file to match other files
allped$sire.pop[which(allped$sire.pop=="ml")]="m"
allped$dam.pop[which(allped$dam.pop=="ml")]="m"#"m"is 5th in the factorized, malinalco in the coancestry matrix is row/column 5. "ML" would come out 6th.
allped$dam.pop<-as.numeric(as.factor(allped$dam.pop))
allped$sire.pop<- as.numeric(as.factor(allped$sire.pop))
#now ped pop numbers reflect alphabet, although dams still reflect field sample order, and sire #s are unrelated and unique for all
# seeds in the same biota are extremely unlikely to share a father, since they come from different infructescences
# note that between biota, seeds sharing a mom may also share a father, since they often come from the same set of infructescences in each biota, but this is not reflected in the pedigree
#all analyses downstream that rely on relatedness are conducted separately across biota treatments for this reason.

envdat <- Apopenvdat[allped$dam.pop,]
cenvdat <- scale(envdat)
#since dam.pop now reflects alphabetical order, can pull out env info this way

# ####MODELS

alldat <- cbind(cpheno,allped,soiltrt,cenvdat)
alldat$dam <- as.character(alldat$dam)
alldat$dam.pop <- as.character(alldat$dam.pop)
traitcols <- c("flowerday","germday","espigalength","shoot","root","finalheight","finalstemw","width3leaf","green")#functional names, matching phenotype data frame

##INTERCEPT
clusterlist <- fitmodelset(traitcols,"~1",alldat)
save(clusterlist,file="~/GenoPheno/traitgxemodelsINT.Rdata")
# ELEVATION
clusterlist <- fitmodelset(traitcols,"~ elevation",alldat)
save(clusterlist,file="~/GenoPheno/traitgxemodelsELEV.Rdata")
#MAT
clusterlist <- fitmodelset(traitcols,"~ TAnn",alldat)
save(clusterlist,file="~/GenoPheno/traitgxemodelsTAnn.Rdata")
# Pann
clusterlist <- fitmodelset(traitcols,"~ Pann",alldat)
save(clusterlist,file="~/GenoPheno/traitgxemodelsPann.Rdata")
#SWC 
clusterlist <- fitmodelset(traitcols,"~ CapacidadDeCampoPercent",alldat)
save(clusterlist,file="~/GenoPheno/traitgxemodelsSWC.Rdata")

traitcols <- c("flowerday","germday","espigalength","shoot","root","finalheight","finalstemw","width3leaf","green")#functional names, matching phenotype data frame
load("~/GenoPheno/traitgxemodelsINT.Rdata")
diccl <- c()
counter <- 1
for(k in 1:9){#traits inside lists of model type
for(j in 1:5){#model type is the first level of list
for(i in 1:10){#repeats innermost level
	diccl[counter] <- clusterlist[[j]][[k]][[i]]$DIC
	counter <- counter +1
	}
}
}
rm(clusterlist)
arraydicclI<-array(diccl,dim=c(10,5,9))#fills first array first(last dimension), first row first(first dimension)

MATmods <- mod.summs("~/GenoPheno/traitgxemodelsTAnn.Rdata",9,5,10)
TAPmods <- mod.summs("~/GenoPheno/traitgxemodelsPann.Rdata",9,5,10)
SWCmods <- mod.summs("~/GenoPheno/traitgxemodelsSWC.Rdata",9,5,10)
ELVmods <- mod.summs("~/GenoPheno/traitgxemodelsELEV.Rdata",9,5,10)

traitcols <- c("flowerday","germday","espigalength","shoot","root","finalheight","finalstemw","width3leaf","green")#functional names, matching phenotype data frame


# ##SUMMARY
randtypes <- c("fp","fps","fpspxs","fpsfxs","fpspxsfxs")#same order as listed in for the cluster list

dicdf <- rbind( data.frame(apply(arraydicclI,c(2,3),mean)), data.frame(apply(MATmods[[1]],c(2,3),mean)), data.frame(apply(TAPmods[[1]],c(2,3),mean)),
		data.frame(apply(SWCmods[[1]],c(2,3),mean)), data.frame(apply(ELVmods[[1]],c(2,3),mean)) )
colnames(dicdf) <- paste(traitcols,"DIC",sep=".")
pvaldf <- rbind(data.frame(matrix(NA,ncol=9,nrow=5)), data.frame(apply(MATmods[[2]],c(2,3),mean)), data.frame(apply(TAPmods[[2]],c(2,3),mean)),
	  data.frame(apply(SWCmods[[2]],c(2,3),mean)), data.frame(apply(ELVmods[[2]],c(2,3),mean)) )
colnames(pvaldf) <- paste(traitcols,"pval",sep=".")
slopedf <- rbind(data.frame(matrix(NA,ncol=9,nrow=5)), data.frame(apply(MATmods[[3]],c(2,3),mean)), data.frame(apply(TAPmods[[3]],c(2,3),mean)),
	  data.frame(apply(SWCmods[[3]],c(2,3),mean)), data.frame(apply(ELVmods[[3]],c(2,3),mean)) )
colnames(slopedf) <- paste(traitcols,"slp",sep=".")

resdf <- data.frame(cbind(dicdf,pvaldf,slopedf))
resdf$envvar <- rep(c("I","MAT","TAP","SWC","ELV"), each = nrow(resdf)/5 )
resdf$randtype <- rep(randtypes, times = 5 )
write.csv(resdf,"~/GenoPheno/04_02_gxeresults.csv",row.names=F)

sink("~/GenoPheno/04_02_gxeresults.txt")
for(i in 1:9){
	print(head(resdf[order(resdf[,i]),c(i,i+9,i+18,ncol(resdf)-1,ncol(resdf))]))
}
sink()
