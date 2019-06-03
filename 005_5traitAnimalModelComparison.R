library(MCMCglmm)

soils <- read.csv("~/GenoPheno/covmat.csv",stringsAsFactors=F,header=T)

allpheno <- read.csv("~/GenoPheno/phenomat.csv",stringsAsFactors=F,header=T)
allpheno <- allpheno[,-c(3,10:12,14,16)]#cuts problematic traits
greens <- (allpheno$green)/(allpheno$red + allpheno$blue + allpheno$green)
allpheno$green <- greens

colnames(allpheno)[1]<- "ID"
#split data into three matrices by soil
traits.mt <- allpheno[which(soils$treatment==3),c(1,2,3,5,6,8)]
traits.mc <- allpheno[which(soils$treatment==5),c(1,2,3,5,6,8)]
traits.tc <- allpheno[which(soils$treatment==7),c(1,2,3,5,6,8)]

allped <- read.csv("~/GenoPheno/pedmat.csv",stringsAsFactors=F,header=T)
#create or read pedigree. no # sharing between dam and sire.#300 plants. 10 moms per pop, 3 sibs per mom, and organized perfectly in trait data, this is what you need
allped$sire.pop[which(allped$sire.pop=="ml")]="m"
allped$dam.pop[which(allped$dam.pop=="ml")]="m"#"m"is 5th.  "ML" would come out 6th. this is for consistency with other scripts/files.
allped$dam.pop<-as.numeric(as.factor(allped$dam.pop))
allped$sire.pop<- as.numeric(as.factor(allped$sire.pop))
#now ped pop numbers reflect alphabet, although dams still reflect field sample order, and sire #s are unrelated and unique for all
# seeds in the same biota are extremely unlikely to share a father, since they come from different infructescences
# note that between biota, seeds sharing a mom may also share a father, since they often come from the same set of infructescences in each biota, but this is not reflected in the pedigree
#all analyses that rely on relatedness are conducted separately across biota treatments for this reason.

#quick function to scale a dataframe that also has an identifier column that should not be scaled
IDscale <- function(inputdf,IDcol=1) {
tmp <- scale(inputdf, center = TRUE, scale = TRUE)#
tmp[,IDcol]<- inputdf[,IDcol]#restore id
return(tmp)
}


##############G MATRICES 
# 
# ###############BY POPxSOIL################# 
#subset allpheno and allped into popxsoil specific subsets
psi <- cbind(rep(1:10,times=3),rep(c(3,5,7),each=10)) #soil trts 3 = mt 5=mc 7 = tc. soils are not alphabetical by numbers, but pop numbers ARE alphabetical
popsoil.pheno <- lapply(1:nrow(psi), function(z) allpheno[which(allped$dam.pop==psi[z,1] & soils$treatment==psi[z,2]),])
popsoil.ped <- lapply(1:nrow(psi), function(z) allped[which(allped$dam.pop==psi[z,1] & soils$treatment==psi[z,2]),])
#scale, center
popsoil.phenoC <-  lapply(1:30,function(z) as.data.frame(IDscale(popsoil.pheno[[z]])))

#format
glm.popsoil.ped <- lapply(1:30, function(z) as.matrix(data.frame(
				id=c(unique(popsoil.ped[[z]]$dam), popsoil.ped[[z]]$sire,popsoil.ped[[z]]$id),
				sire = c(rep(NA,length.out=40),popsoil.ped[[z]]$sire), dam=c(rep(NA,length.out=40),popsoil.ped[[z]]$dam)))
				)#no info on dams or sires, only offspring. here, 30 offspring, 10 moms, 30 dads
for(i in 1:30){colnames(glm.popsoil.ped[[i]])[1] <- "animal"}

glm.popsoil.pheno  <- popsoil.phenoC
 for(i in 1:30){
 	colnames(glm.popsoil.pheno[[i]])[1] <- "animal" } #dampop column unnecessary, only one pop!	

#prior
priorpopsoil<-list(G=list(G1=list(V=diag(5)*0.02,nu=6)), R=list(V=diag(5)*0.02,nu=6))#list(G1= ..., G2=list(V=diag(9)*0.02,nu=10)) #(this exactly the same as priorpop)
priorpop<-list(G=list(G1=list(V=diag(5)*0.02,nu=6)), R=list(V=diag(5)*0.02,nu=6))#V = dim(V)*.02 nu = dim(V) + 1 
priordownV<-list(G=list(G1=list(V=diag(5)*0.01,nu=6)), R=list(V=diag(5)*0.01,nu=6))#approaching V = dim(V)*0 nu = dim(V) +1 ##cant run this or below bc not posdef.
priordownnu<-list(G=list(G1=list(V=diag(5)*0.01,nu=2)), R=list(V=diag(5)*0.01,nu=2))#approaching V = dim(V)*0 nu = dim(V) - 3
#all three of the above priors are from Hadfield course notes "3.6.1 Priors for us structures" and is a hopefully weakly informative, proper prior of the inverse wishart.
#Hadfield source is cited in main text.

save(glm.popsoil.ped,file="~/GenoPheno/Popsoil_Peds.Rdata")

#models
#CHOSEN PROPER weakly informative prior.
popsoil5.animal.mods <- lapply(1:30, function(z) MCMCglmm(cbind(flowerday,germday,shoot,root,finalstemw) 
~ 1 , random=~us(trait):animal, rcov=~us(trait):units, #note removal of trait main effect, should be 0 since data is centered, appears to slow fitting somewhat
	family=c("gaussian","gaussian","gaussian","gaussian","gaussian"),
	pedigree=glm.popsoil.ped[[z]],data=glm.popsoil.pheno[[z]],prior=priorpopsoil,verbose=FALSE,pr=TRUE,nitt=1000000, thin=1000, burnin=100000))#
save(popsoil5.animal.mods,file="~/GenoPheno/ComparisonG/popsoil.animal.mods.5trait.Rdata")
##nearly improper and nearly non-informative priors.
popsoil5.animal.mods.imp <- lapply(1:30, function(z) MCMCglmm(cbind(flowerday,germday,shoot,root,finalstemw) 
~ 1 , random=~us(trait):animal, rcov=~us(trait):units, #note removal of trait main effect, should be 0 since data is centered, appears to slow fitting somewhat
	family=c("gaussian","gaussian","gaussian","gaussian","gaussian"),
	pedigree=glm.popsoil.ped[[z]],data=glm.popsoil.pheno[[z]],prior=priordownV,verbose=FALSE,pr=TRUE,nitt=1000000, thin=1000, burnin=100000))#
save(popsoil5.animal.mods.imp,file="~/GenoPheno/ComparisonG/popsoil.animal.mods.nearlyimp.prior.5trait.Rdata")
popsoil5.animal.mods.ninf <- lapply(1:30, function(z) MCMCglmm(cbind(flowerday,germday,shoot,root,finalstemw) 
~ 1 , random=~us(trait):animal, rcov=~us(trait):units, #note removal of trait main effect, should be 0 since data is centered, appears to slow fitting somewhat
	family=c("gaussian","gaussian","gaussian","gaussian","gaussian"),
	pedigree=glm.popsoil.ped[[z]],data=glm.popsoil.pheno[[z]],prior=priordownnu,verbose=FALSE,pr=TRUE,nitt=1000000, thin=1000, burnin=100000))#
save(popsoil5.animal.mods.ninf,file="~/GenoPheno/ComparisonG/popsoil.animal.mods.nearlynoninf.prior.5trait.Rdata")
