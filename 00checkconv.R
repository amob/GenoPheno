
library(MCMCglmm)#dependencies

make.array <- function(model.list,ntraits,nsamples,whichtype){ #whichtype , for DS this can be done across soils only, I believe.
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

traceallG <- function(GmatArray,outfile,ntraits,nmats,traitnames){
for(k in 1:nmats){
	pdf(paste(outfile,k,".pdf",sep=""),width=ntraits*3,height=ntraits*3)
	par(mfrow=c(ntraits,ntraits),mar=c(4,6,5,3))
	for(i in 1:ntraits){
		for(j in 1:ntraits){
		plot(GmatArray[i,j,k,],cex.lab=2.4,cex.axis=2.4,cex.main=2.4,ylab=traitnames[i],main=traitnames[j])
		}
	
	}
	dev.off()
}}

tnames5 <- c("Days to flowering","Days to germination","Shoot biomass","Root biomass","Stem width")
load("~/GenoPheno/ComparisonG/popsoil.animal.mods.5trait.Rdata")
An.popsoil <- make.array( popsoil5.animal.mods , 5 , nrow((popsoil5.animal.mods[[1]])$VCV) , whichtype="animal" )
traceallG(An.popsoil,"~/GenoPheno/AnimalConvergenceChecks/properprior5traits/TESTprop5AnimalPopSoiltraitCheckConvcombo_",5,30,tnames5)
rm(popsoil5.animal.mods);rm(An.popsoil)

load("~/GenoPheno/ComparisonG/popsoil.animal.mods.nearlyimp.prior.5trait.Rdata")
An.popsoil <- make.array( popsoil5.animal.mods.imp , 5 , nrow((popsoil5.animal.mods.imp[[1]])$VCV) , whichtype="animal" )
traceallG(An.popsoil,"~/GenoPheno/AnimalConvergenceChecks/nearlyimproperprior5traits/improp5AnimalPopSoiltraitCheckConvcombo_",5,30)
rm(popsoil5.animal.mods.imp);rm(An.popsoil)

load("~/GenoPheno/ComparisonG/popsoil.animal.mods.nearlynoninf.prior.5trait.Rdata")
An.popsoil <- make.array( popsoil5.animal.mods.ninf , 5 , nrow((popsoil5.animal.mods.ninf[[1]])$VCV) , whichtype="animal" )
traceallG(An.popsoil,"~/GenoPheno/AnimalConvergenceChecks/nearlynoninfprior5traits/noninf5AnimalPopSoiltraitCheckConvcombo_",5,30)
rm(popsoil5.animal.mods.ninf);rm(An.popsoil)
