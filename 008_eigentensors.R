#purpose: eigentensor analysis of g matrices for popsoil, pop, soil, sets of matrices. compared to randomized matrices

library(gdata)
library(matrixcalc)
library(MCMCglmm)#dependencies

range01=function(x){
newnums=(x-min(x,na.rm=T))/(max(x,na.rm=T)-min(x,na.rm=T))
return(newnums)
}

############this function adapted from  Aguirre et al 2014
#output arranged in [trait,trait,pop,mcmc] fashion
make.array <- function(model.list,ntraits,nsamples,whichtype){ #whichtype, allows using on ancestral G matrices estimated in driftsel, we do not report these results
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

####MAKE RANDOM ARRAYS.
#######Random array code, adapted from Aguirre et al 2014, modifications include making it a function,looping over pop/popsoil/soil G matrices
#requires pedigrees with unique codes for individuals that are not shared across Gs
#if popefs is TRUE, then popmeans must be a list of dataframes with fixed effects groups in columns and samples from model output in rows
make.rand <- function(model.array, ped.list,n,m,popefs=FALSE,popmeans=NULL,popid.list=NULL,individs=NULL){#individs now req'd, and is # individuals used to estimate each G(balanced)
	MCMCsamp <- dim(model.array)[4]
	rand.Garray <- array(,c(n,n,m,MCMCsamp))
	a.pop <- cbind(seq(from=1,by=individs,length.out=length(ped.list)),seq(from=individs,by=individs,length.out=length(ped.list)))
	for (i in 1:MCMCsamp){
    	if(popefs==TRUE) { pops <- popid.list[[1]]} else {pops <- NULL}
    	if(popefs==TRUE) {groupefs <- popmeans[[1]][,,i]} else{groupefs <-  NULL}
    	rbv.nodes <- !is.na(ped.list[[1]][,2])
    	pop.bv <- rbv(ped.list[[1]],model.array[,,1,i],ggroups= pops, gmeans= groupefs)[rbv.nodes,] #rbv function from MCMCglmm. change from Aguirre is ADDED POP EFFS
    	for(j in 2:m)  {  
    		if(popefs==TRUE) {pops <- popid.list[[j]]} else {pops <- NULL}
    		if(popefs==TRUE) {groupefs <- popmeans[[j]][,,i]} else{groupefs <-  NULL}
    		rbv.nodes <- !is.na(ped.list[[j]][,2])
    		pop.bv <- rbind( pop.bv, rbv(ped.list[[j]],model.array[,,j,i],ggroups=pops,gmeans=groupefs)[rbv.nodes,] )  }#
		for(j in 1:m) {
		rand.pop.bv <- pop.bv[sample(dim(pop.bv)[1],replace=F),]
		rand.Garray[,,j,i] <- cov(rand.pop.bv[a.pop[j,1]:a.pop[j,2],])#
		}#
  	}
  return(rand.Garray)
}

######tensor analysis, Aguierre et al 2014 code. AMOB edited this only in defining input variables in-function.
#set function
covtensor <- function(Gs){
    if (dim(Gs)[[1]] != dim(Gs)[[2]]){
      stop("G array must be of order n x n x m x MCMCsamp")
    }
    if (is.na(dim(Gs)[4])) {
      stop("There are no MCMCsamples")
    }
    m <- dim(Gs)[3]##AMOB edits set m and n and MCMCsamp variables
    n <- dim(Gs)[1]##AMOB, ditto
    MCMCsamp <- dim(Gs)[4]##AMOB, ditto    
    neigten <- n*(n+1)/2 
    #Number of eigentensors
    MCMC.S <- array(,c(neigten, neigten, MCMCsamp))
    dimnames(MCMC.S) <- list(paste("e", 1:neigten, sep=""), paste("e", 1:neigten, sep=""))
      for (k in 1:MCMCsamp){
        MCMCG <- Gs[,,,k] 
          MCMCvarmat <- t(apply(MCMCG, 3, diag)) 
          #find the variances of the kth G and store them 
          MCMCcovmat <- t(apply(MCMCG, 3, lowerTriangle)) 
          #find the covariances of the kth G and store them
          MCMC.S[1:n,1:n, k] <- cov(MCMCvarmat, MCMCvarmat) 
          #fill the upper left quadrant of the kth S
          MCMC.S[(n+1):neigten,(n+1):neigten, k] <- 2*cov(MCMCcovmat, MCMCcovmat)
          #fill the lower right quadrant of the kth S
          MCMC.S[1:n,(n+1):neigten, k] <- sqrt(2)*cov(MCMCvarmat, MCMCcovmat)
          #fill the upper right quadrant of the kth S
          MCMC.S[(n+1):neigten,1:n, k] <- sqrt(2)*cov(MCMCcovmat, MCMCvarmat)
          #fill the lower left quadrant of the kthS
        }  
        av.S <- apply(MCMC.S, 1:2, mean)
        #posterior mean S
        av.S.val <- eigen(av.S)$values
        #eigenvalues of posterior mean S 
        av.S.vec <- eigen(av.S)$vectors
        #eigenvalues of posterior mean S
        eTmat <- array(, c(n, n, neigten))
        dimnames(eTmat) <- list(traitnames, traitnames, paste("E", 1:neigten, sep=""))  
          for (i in 1:neigten){
            emat <- matrix(0, n, n) 
            lowerTriangle(emat) <- 1/sqrt(2)*av.S.vec[(n+1):neigten,i]
            emat <- emat + t(emat)
            diag(emat) <- av.S.vec[1:n,i]
            eTmat[,,i] <- emat 
          }
          #construct the second-order eigentensors of posterior mean S
          eT.eigen <- array(, c(n+1, n, neigten))
            for (i in 1:neigten){
              eT.eigen[1,,i] <- t(eigen(eTmat[,,i])$values) 
              #Eigenvalues of the ith eigentensor
              eT.eigen[2:(n+1),,i] <- eigen(eTmat[,,i])$vectors 
              #Eigenvectors of the ith eigentensor
              eT.eigen[,,i] <- eT.eigen[,order(abs(eT.eigen[1,,i]), decreasing = T), i]
            }
            MCMC.S.val <- matrix(, MCMCsamp, neigten)
            colnames(MCMC.S.val) <- paste("E", 1:neigten, sep="")
            for (i in 1:MCMCsamp){
              for(j in 1:neigten){
                MCMC.S.val[i,j] <- t(av.S.vec[,j]) %*% MCMC.S[,,i] %*% av.S.vec[,j]
              }
            }
            #posterior distribution of the genetic variance for the eigenvectors of posterior mean S
            av.G.coord <- array(, c(m, neigten, 1))
            dimnames(av.G.coord) <- list(Gnames, paste("E", 1:neigten, sep=""))
              for (i in 1:neigten){
                av.G.coord[,i,] <- apply((apply(Gs, 1:3, mean)) , 3, frobenius.prod, y = eTmat[,,i])
              }
              #Coordinates of the jth avG for the eigentensors of posterior mean S
              MCMC.G.coord <- array(, c(m, neigten, MCMCsamp))
              dimnames(MCMC.G.coord) <- list(Gnames, paste("E", 1:neigten, sep=""))
                for (i in 1:neigten){
                  MCMC.G.coord[,i,] <- apply(Gs, 3:4, frobenius.prod, y = eTmat[,,i])
                }
                #Coordinates of the kth MCMC sample of the jth G for the eigentensors of posterior mean S
        tensor.summary <- data.frame(rep(av.S.val,each=n), t(data.frame(eT.eigen)))
        colnames(tensor.summary) <- c("S.eigval", "eT.val", traitnames)
        rownames(tensor.summary)<- paste(paste("e", rep(1:neigten, each=n), sep=""), rep(1:n,neigten), sep=".")
  list(tensor.summary = tensor.summary, av.S = av.S, eTmat = eTmat, av.G.coord = av.G.coord, MCMC.S = MCMC.S, MCMC.S.val = MCMC.S.val, MCMC.G.coord = MCMC.G.coord)
}


######
#read in models
load("~/GenoPheno/ComparisonG/popsoil.animal.mods.5trait.Rdata")
An.popsoil.5T <- make.array( popsoil5.animal.mods , 5 , nrow((popsoil5.animal.mods[[1]])$VCV) , whichtype="animal" )
load("~/GenoPheno/ComparisonG/popsoil.animal.mods.nearlyimp.prior.5trait.Rdata")
An.popsoil.I5T <- make.array( popsoil5.animal.mods.imp , 5 , nrow((popsoil5.animal.mods.imp[[1]])$VCV) , whichtype="animal" )
load("~/GenoPheno/ComparisonG/popsoil.animal.mods.nearlynoninf.prior.5trait.Rdata")
An.popsoil.N5T <- make.array( popsoil5.animal.mods.ninf , 5 , nrow((popsoil5.animal.mods.ninf[[1]])$VCV) , whichtype="animal" )

#read in pedigrees, dataframes in list
load("~/GenoPheno/Popsoil_Peds.Rdata")#glm.popsoil.ped

envmatall <-read.csv("~/GenoPheno/AlphabeticalPopEnvDat.csv")
 
####RUN TENSOR ANALYSIS

#popsoil animal mods were iterated over 3, 5, 7 sample order of soil sites. but alphabetical order of plant pops
randAn.popsoil.5T <- make.rand(An.popsoil.5T,	glm.popsoil.ped,5,30,individs=30)#5 traits, 30 pops, 30 individuals per pop
randAn.popsoil.I5T <- make.rand(An.popsoil.I5T,	glm.popsoil.ped,5,30,individs=30)
randAn.popsoil.N5T <- make.rand(An.popsoil.N5T,	glm.popsoil.ped,5,30,individs=30)

tnames5 <- c("Days to flowering","Days to germination","Shoot biomass","Root biomass","Stem width")#covtensor calls for traitnamesvector
Gnames.ps <- paste(rep(c("cl","cu","da","fp","m","mc","mt","tc","tx","tz"),times=3),".",rep(c("mc","mt","tc"),each=10),sep="")
Gnames <- Gnames.ps #covtensor requires Gnames and traitnames objects to exist
traitnames <- tnames5
tens.PS.5T <- covtensor(An.popsoil.5T)
tens.PS.I5T <- covtensor(An.popsoil.I5T)
tens.PS.N5T <- covtensor(An.popsoil.N5T)
#has 15 tensors bc 15 gmat elements, which is less than 30 matrices

randtens.PS.5T <- covtensor(randAn.popsoil.5T)
randtens.PS.I5T <- covtensor(randAn.popsoil.I5T)
randtens.PS.N5T <- covtensor(randAn.popsoil.N5T)


##HPDeT function below also adapted from Aguirre code:  HPDIs around eigenvalues of tensors
 HPDeT <- function(tensorS,randtensorS) {
 			return( round(cbind(HPDinterval(as.mcmc(tensorS), prob=0.95), 
 					HPDinterval(as.mcmc(randtensorS), prob=0.95)),digits=3) ) 
 			}
HPDet5T<-HPDeT(as.mcmc(tens.PS.5T$MCMC.S.val),as.mcmc(randtens.PS.5T$MCMC.S.val))
HPDet5T# tensors 13-15 NS. first two fairly sig.
HPDeT(as.mcmc(tens.PS.I5T$MCMC.S.val),as.mcmc(randtens.PS.I5T$MCMC.S.val))#only tensors 12,13 & 15 ns
HPDeT(as.mcmc(tens.PS.N5T$MCMC.S.val),as.mcmc(randtens.PS.N5T$MCMC.S.val))# tensors 1, 9, 12, 13, 15.
#eigentensor 1 has more error from G matrices fit with nearly noninformative prior

pdf("~/GenoPheno/gmatrix/EigenValS_AnPopSoilG5trait.pdf",height=4,width=4)
plot(unique(tens.PS.5T$tensor.summary[,1])~c(1:15),pch=16,ylab="Eigenvalues of eigentensors",xlab="",xaxt="n",main="Populations",ylim=c(0,2.5),xlim=c(0.5,16),cex=1,cex.axis=1,cex.lab=1,cex.main=1.2)
points(unique(randtens.PS.5T$tensor.summary[,1])~c(1:15 + .5),pch=1,cex=1)
arrows(c(1:15),unique(tens.PS.5T$tensor.summary[,1])[1:15],x1=c(1:15),y1=HPDet5T[1:15,1],length = 0)
arrows(c(1:15),unique(tens.PS.5T$tensor.summary[,1])[1:15],x1=c(1:15),y1=HPDet5T[1:15,2],length = 0)
arrows(1:15 + .5,unique(randtens.PS.5T$tensor.summary[,1])[1:15],x1=c(1:15 + .5),y1=HPDet5T[1:15,3],length = 0)
arrows(c(1:15 + .5),unique(randtens.PS.5T$tensor.summary[,1])[1:15],x1=c(1:15 + .5),y1=HPDet5T[1:15,4],length = 0)
dev.off()

#variance explained
varexplnd <- function(tensorsum, traitnum){
	tensvar <- unique(tensorsum[,1])/sum(unique(tensorsum[,1]))
	tens2vec2 <- sapply(seq(from=1,by=traitnum,length.out=2), function(z) abs(tensorsum[z:(z+traitnum-1),2]/sum(abs(tensorsum[z:(z+traitnum-1),2]))) ) 
	return(list(tensvar,tens2vec2))
	}
varexplnd(tens.PS.5T$tensor.summary,5)
varexplnd(tens.PS.I5T$tensor.summary,5)
varexplnd(tens.PS.N5T$tensor.summary,5)

rbd<-colorRampPalette(c(rgb(.8,0,0),rgb(.8,0,0,alpha=.5),rgb(1,1,1),rgb(0,0,.8,alpha=.5),rgb(0,0,.8)))

pdf("~/GenoPheno/gmatrix/TraitsPopSoilEigenVofEigenT5traits.pdf",height=3,width=3)
par(mar=c(8,3,1,5))
image(1:5,c(.5,1.5,2.5,3.5),t(as.matrix(tens.PS.5T$tensor.summary[c(1,2,6,7),3:7])),ylim=c(0,4),zlim=c(-1,1),axes=F,ylab="",xlab="",col=rbd(9))
points(0,0)# so that labels will display properly
axis(side=1, at=1:5, labels=tnames5, las=2, cex.axis = 0.8)
axis(side=2, at=c(.5,1.5,2.5,3.5), labels=c("T1V1","T1V2","T2V1","T2V2"), las=2, cex.axis = 0.8)
dev.off()

tens.PS.5T$eTmat[,,1] - tens.PS.I5T$eTmat[,,1] # very small changes, see below
tens.PS.5T$eTmat[,,1] - tens.PS.N5T$eTmat[,,1] #changes 0.01 - 0.02;; most well below .01

ltens1ps5 <- tens.PS.5T$eTmat[,,1]
ltens1ps5[lower.tri(ltens1ps5)] <- 0
pdf("~/GenoPheno/gmatrix/PopsoilTensor1.5trait.pdf",height=4,width=4)
par(mar=c(9,9,1,1))
image(1:5,seq(from = .5, to =5.5, by =1),((ltens1ps5)[,5:1]),ylim=c(0.5,5.5),axes=F,ylab="",xlab="",zlim=c(-.5,.5),col=rbd(151))
axis(side=1, at=1:5, labels=(tnames5), las=2, cex.axis = 1)
axis(side=2, at=1:5, labels=rev(tnames5), las=2, cex.axis = 1)
dev.off()
