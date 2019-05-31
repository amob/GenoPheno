
library(MCMCglmm)#dependencies

############code/structures adapted from  Aguirre et al 2014
#to do this I need output arranged in [trait,trait,pop,mcmc] fashion
make.array <- function(model.list,ntraits,nsamples,whichtype){ #whichtype allows input of a driftsel object, not currently reported
	mcmc.array <- array(,dim=c(ntraits,ntraits,length(model.list),nsamples))
	if(whichtype=="DS"){
		for(i in 1:length(model.list)){
			mcmc.array[,,i,] <- model.list[[i]]$G
			}
	} else if(whichtype=="animal"){
		gcols <- grep(colnames(model.list[[1]]$VCV),pattern="animal")
		for(i in 1:length(model.list)){
			for(j in 1:nsamples){
				mcmc.array[,,i,j] <- matrix(model.list[[i]]$VCV[j,gcols],ncol=ntraits)
			}
		}
	} else {print("whichtype should equal \"DS\" or \"Animal\"")}
	return(mcmc.array)
}
###########
######
#read in models
load("~/GenoPheno/ComparisonG/popsoil.animal.mods.5trait.Rdata")
An.popsoil.5T <- make.array( popsoil5.animal.mods , 5 , nrow((popsoil5.animal.mods[[1]])$VCV) , whichtype="animal" )
load("~/GenoPheno/ComparisonG/popsoil.animal.mods.5traitSQ.Rdata")
An.popsoil.5TSQ <- make.array( popsoil5.animal.modsSQ , 5 , nrow((popsoil5.animal.modsSQ[[1]])$VCV) , whichtype="animal" )
load("~/GenoPheno/ComparisonG/popsoil.animal.mods.nearlyimp.prior.5trait.Rdata")
An.popsoil.I5T <- make.array( popsoil5.animal.mods.imp , 5 , nrow((popsoil5.animal.mods.imp[[1]])$VCV) , whichtype="animal" )

envmatall <-read.csv("~/GenoPheno/AlphabeticalPopEnvDat.csv")

traitnames <- c("DaystoFlower","DaystoGermination","TasselLength","ShootMass","RootMass","Height","StemWidth","LeafWidth","StemColor")#covtensor calls for 
tnames5 <- c("Days to flowering","Days to germination","Shoot biomass","Root biomass","Stem width")#covtensor calls for traitnamesvector


rbd<-colorRampPalette(c(rgb(.8,0,0),rgb(.8,0,0,alpha=.5),rgb(1,1,1),rgb(0,0,.8,alpha=.5),rgb(0,0,.8)))

Gnames.ps <- paste(rep(c("cl","cu","da","fp","m","mc","mt","tc","tx","tz"),times=3),".",rep(c("mt","mc","tc"),each=10),sep="")#represents actual order of Gs. mt mc swapped here relative to some other arrangments.
#Gnames.ps useful for checking rearranging vectors to make sure they rearrange things as you think
swapsoilpoptemp <- rep(c(10,5,9,7,3,4,6,1,8,2),times=3) + rep(c(0,10,20),each=10)#its going to fill by rows below, so we want first one soil, then the next then the last, almost as it is., pops in temp order
AnPopsoilGs.5T <- apply(An.popsoil.5T,c(1,2,3),mean) [,,swapsoilpoptemp]#all three tz matrices first, all three cu last, rearranged to be mc mt tc alphabetical
for(i in 1:30){AnPopsoilGs.5T[,,i] <- (AnPopsoilGs.5T[,,i])*lower.tri(AnPopsoilGs.5T[,,i],diag=T)}


pdf("~/GenoPheno/gmatrix/RawGAnPopsoil5trait.pdf",height=7,width=24)
par(mfrow=c(3,10))#
par(mar=c(1,1,1,2))
for(i in 1:30){
image(1:5,seq(from=.5,by=1, length.out=5),apply(AnPopsoilGs.5T[,,i],2,rev),ylim=c(0,6),zlim=c(-1,1),axes=F,ylab="",xlab="",col=rbd(101))
}
dev.off()
