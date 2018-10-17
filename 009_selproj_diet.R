#purpose: selection skewers analysis of g matrices
#library(gdata)
#library(matrixcalc)
library(MCMCglmm)#dependencies
#library(LaplacesDemon)
#library(ellipse)
#library(MASS)
library(ks)#installed with dependencies

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

PhenSelection <- function(Gs,phensel){# phensel = c(-1,0,0,0,0,0,0,0,0)
 	n <- dim(Gs)[[1]]
 	m <- dim(Gs)[[3]]
 	MCMCsamp <- dim(Gs)[[4]]
 	 phen.vec <- phensel/(sqrt(sum(phensel^2)))
		deltaZs <- array(,c(MCMCsamp,n, m))
		   for(k in 1:m){
    	      deltaZs[,,k]<- t(sapply(1:MCMCsamp, function(z) phen.vec%*%Gs[,,k,z]))
        	}
    	return(deltaZs)
}


load("~/GenoPheno/ComparisonG/popsoil.animal.mods.5trait.Rdata")
An.popsoil.5T <- make.array( popsoil5.animal.mods , 5 , nrow((popsoil5.animal.mods[[1]])$VCV) , whichtype="animal" )
load("~/GenoPheno/ComparisonG/popsoil.animal.mods.nearlyimp.prior.5trait.Rdata")
An.popsoil.I5T <- make.array( popsoil5.animal.mods.imp , 5 , nrow((popsoil5.animal.mods.imp[[1]])$VCV) , whichtype="animal" )
load("~/GenoPheno/ComparisonG/popsoil.animal.mods.nearlynoninf.prior.5trait.Rdata")
An.popsoil.N5T <- make.array( popsoil5.animal.mods.ninf , 5 , nrow((popsoil5.animal.mods.ninf[[1]])$VCV) , whichtype="animal" )


#read in pedigrees, dataframes in list
load("~/GenoPheno/Popsoil_Peds.Rdata")#glm.popsoil.ped

envmatall <-read.csv("~/GenoPheno/AlphabeticalPopEnvDat.csv")

tnames5 <- c("Days to flowering","Days to germination","Shoot biomass","Root biomass","Stem width")#covtensor calls for traitnamesvector


##popsoil
#popsoil animal mods were iterated over 3, 5, 7 sample order of soil sites. but alphabetical order of plant pops
Gnames.ps <- paste(rep(c("cu","cl","da","fp","m","mc","mt","tc","tx","tz"),times=3),".",rep(c("mc","mt","tc"),each=10),sep="")#alphabetical.

#skewers
#sel.Anpopsoil <- PhenSelection(An.popsoil,c(1,0,0,0,0,0,0,0,0))
sel.Anpopsoil.5T <- PhenSelection(An.popsoil.5T,c(1,0,0,0,0))
sel.Anpopsoil.I5T <- PhenSelection(An.popsoil.I5T,c(1,0,0,0,0))
sel.Anpopsoil.N5T <- PhenSelection(An.popsoil.N5T,c(1,0,0,0,0))#these are now only deltaZs

#plots
rbd<-colorRampPalette(c(rgb(.8,0,0),rgb(.8,0,0,alpha=.5),rgb(1,1,1),rgb(0,0,0.8,alpha=.5),rgb(0,0,0.8)))

Gnames.ps <- paste(rep(c("cu","cl","da","fp","m","mc","mt","tc","tx","tz"),times=3),".",rep(c("mt","mc","tc"),each=10),sep="")#represents actual order of Gs. 
#Gnames.ps useful for checking rearranging vectors to make sure they rearrange things as you think
swapsoilpoptemp <- rep(c(10,5,9,7,3,4,6,1,8,2),times=3) + rep(c(0,10,20),each=10)#its going to fill by rows below, so we want first one soil, then the next then the last, almost as it is., pops in temp order
AnPopsoilGs.5T <- apply(An.popsoil.5T,c(1,2,3),mean) [,,swapsoilpoptemp]#all three tz matrices first, all three cu last, rearranged to be soil treatment mc mt tc alphabetical
for(i in 1:30){AnPopsoilGs.5T[,,i] <- (AnPopsoilGs.5T[,,i])*lower.tri(AnPopsoilGs.5T[,,i],diag=T)}

contcols1 <- rgb(0,0,0,alpha=seq(from = 0,to=1,length.out=5))
contcols2 <- rgb(0,1,0,alpha=seq(from = 0,to=1,length.out=5))
bandwithmat <- matrix(c(0.02,.008,.008,0.02),nrow=2)

#when AnPopsoilGs is ordered by swapsoilpoptemp, we want in order of temp for both...mt is the 4th warmest pop, mc the 7th, tc the 9th.
#so, plot 4, 7, 9. first row. 14,17,19, second row, 24,27,29 third row Gnames.ps[swapsoilpoptemp][c(4,7,9,14,17,19,24,27,29)], rows are biota, colunmns are pops
pdf("~/GenoPheno/cont2SelSkews3_5Trait.pdf",height=13,width=12)
layout(matrix(1:16,ncol=4,byrow=T),heights=c(8,8,8,4),widths=c(4,8,8,8))
par(mar=c(3,3,7,3))
plot(NA,NA,bty="n",xlim=c(0,1),ylim=c(0,1),yaxt="n",xaxt="n")
image(1:5,seq(from=.5,by=1, length.out=5),apply(AnPopsoilGs.5T[,,12],2,rev),ylim=c(0,5),zlim=c(-1,1),axes=F,ylab="",xlab="",col=rbd(101))
points(0,0)# now these will display properly
axis(side=2, at=(1:5)-.5, labels=rev(tnames5), las=2, cex.axis = 2)
mtext("Malinalco plants",side=3,line=1,cex=1.5)
mtext("Biota15.0",side=3, line=3,adj=1.55,cex=1.5)
image(1:5,seq(from=.5,by=1, length.out=5),apply(AnPopsoilGs.5T[,,18],2,rev),ylim=c(0,5),zlim=c(-1,1),axes=F,ylab="",xlab="",col=rbd(101))
points(0,0)# now these will display properly
mtext("Lower Calimaya Plants",side=3,line=1,cex=1.5)
kde11 <- kde(sel.Anpopsoil.5T[,c(1,4),11], H=bandwithmat, approx.cont=T)
kde15 <- kde(sel.Anpopsoil.5T[,c(1,4),15], H=bandwithmat, approx.cont=T)
plot(kde11,cont=c(50,80,90,95),col=contcols1[5:2],ylim=c(-1,1.5),xlim=c(-0.2,2),drawpoints=T,pch=16,cex=.5,col.pt=rgb(0,0,0,alpha=.1))
plot(kde15,cont=c(50,80,90,95),col=contcols2[5:2],add=T,drawpoints=T,pch=16,cex=.5,col.pt=rgb(0,1,0,alpha=.1))
mtext("Root Mass Response",side=2,line=3,cex=1.5)
mtext("Flowering Time Response",side=1,line=4,cex=1.5)
mtext("Selection for later flowering",side=3,line=2,cex=1.4)#
legend(-3.75,.75,c("Malinalco","Calimaya"),col=c(rgb(0,1,0),rgb(0,0,0)),pch=16,bty="n")

plot(NA,NA,bty="n",xlim=c(0,1),ylim=c(0,1),yaxt="n",xaxt="n")
image(1:5,seq(from=.5,by=1, length.out=5),apply(AnPopsoilGs.5T[,,2],2,rev),ylim=c(0,5),zlim=c(-1,1),axes=F,ylab="",xlab="",col=rbd(101))
axis(side=2, at=(1:5)-.5, labels=rev(tnames5), las=2, cex.axis = 2)
mtext("Malinalco plants",side=3,line=1,cex=1.5)
mtext("Biota14.3", side = 3,line=3,adj=1.55,cex=1.5)
image(1:5,seq(from=.5,by=1, length.out=5),apply(AnPopsoilGs.5T[,,8],2,rev),ylim=c(0,5),zlim=c(-1,1),axes=F,ylab="",xlab="",col=rbd(101))
mtext("Lower Calimaya Plants",side=3,line=1,cex=1.5)
kde1 <- kde(sel.Anpopsoil.5T[,c(1,4),1],H=bandwithmat)
kde5 <- kde(sel.Anpopsoil.5T[,c(1,4),5],H=bandwithmat)
plot(kde1,cont=c(50,80,90,95),col=contcols1[5:2],ylim=c(-1,1.5),xlim=c(-0.2,2),drawpoints=T,pch=16,cex=.5,col.pt=rgb(0,0,0,alpha=.1))
plot(kde5,cont=c(50,80,90,95),col=contcols2[5:2],add=T,drawpoints=T,pch=16,cex=.5,col.pt=rgb(0,1,0,alpha=.1))
mtext("Root Mass Response",side=2,line=3,cex=1.5)
mtext("Flowering Time Response",side=1,line=4,cex=1.5)

plot(NA,NA,bty="n",xlim=c(0,1),ylim=c(0,1),yaxt="n",xaxt="n")
image(1:5,seq(from=.5,by=1, length.out=5),apply(AnPopsoilGs.5T[,,22],2,rev),ylim=c(0,5),zlim=c(-1,1),axes=F,ylab="",xlab="",col=rbd(101))
points(0,0)# now these will display properly
mtext("Malinalco plants",side=3,line=1,cex=1.5)
mtext("Biota13.0",side=3, line=3,adj=1.55,cex=1.5)
axis(side=1, at=1:5, labels=(tnames5), las=2, cex.axis = 2)
axis(side=2, at=(1:5)-.5, labels=rev(tnames5), las=2, cex.axis = 2)
image(1:5,seq(from=.5,by=1, length.out=5),apply(AnPopsoilGs.5T[,,28],2,rev),ylim=c(0,5),zlim=c(-1,1),axes=F,ylab="",xlab="",col=rbd(101))
points(0,0)# now these will display properly
axis(side=1, at=1:5, labels=(tnames5), las=2, cex.axis = 2)
mtext("Lower Calimaya Plants",side=3,line=1,cex=1.5)
kde21 <- kde(sel.Anpopsoil.5T[,c(1,4),21], H=bandwithmat, approx.cont=T)
kde25 <- kde(sel.Anpopsoil.5T[,c(1,4),25], H=bandwithmat, approx.cont=T)
plot(kde21,cont=c(50,80,90,95),col=contcols1[5:2],ylim=c(-1,1.5),xlim=c(-0.2,2),drawpoints=T,pch=16,cex=.5,col.pt=rgb(0,0,0,alpha=.1))
plot(kde25,cont=c(50,80,90,95),col=contcols2[5:2],add=T,drawpoints=T,pch=16,cex=.5,col.pt=rgb(0,1,0,alpha=.1))
mtext("Root Mass Response",side=2,line=3,cex=1.5)
mtext("Flowering Time Response",side=1,line=4,cex=1.5)
legend(0.5,-0.3,c("Malinalco","Calimaya"),col=c(rgb(0,1,0),rgb(0,0,0)),lty=1,bty="n",cex=2)
dev.off()

contcols1 <- rgb(0,0,0,alpha=seq(from = 0,to=1,length.out=5))
contcols2 <- rgb(0,1,0,alpha=seq(from = 0,to=1,length.out=5))
contcols3 <- rgb(1,0,1,alpha=seq(from = 0,to=1,length.out=5))
bandwithmat <- matrix(c(0.02,.008,.008,0.02),nrow=2)

all.sel.contours <- lapply(1:30, function(z) kde(sel.Anpopsoil.5T[,c(1,4),z],H=bandwithmat,approx.cont=T))
popmats <- matrix( rep(c(10,0,20),times=10) + rep(c(10,5,9,7,3,4,6,1,8,2),each=3), ncol=3,byrow=T) # all.sel.contours is alphabetical and needs to be in MAT order
popnamestemp <- c("Tepoztlán 19.8","Malinalco 18.6","Texcoco 15.3","San Mateo Tezoquipan 15.0","Tenango del Aire 14.7",
					"San Francisco Pedregal 14.4","San Matías Cuijingo 14.3","Calimaya Lower 13.2","Toluca 13.0","Calimaya Upper 12.9") #CL CU DA FP M MC MT TC TX TZ
pdf("~/GenoPheno/TESTsoil_sel_all.pdf",height=15,width=7)
layout(matrix(1:10,ncol=2,byrow=T))
par(oma=c(4,4,1,0))
par(mar=c(3,3,2,1))
for(i in 1:10){
	plot(all.sel.contours[[popmats[i,1]]],cont=c(50,80,90,95),col=contcols1[5:2],      drawpoints=F,pch=16,cex=.5,col.pt=rgb(0,0,1,alpha=.1),
		ylim=c(-1,1.75),xlim=c(-0.2,1.75),xlab="",ylab="")
	plot(all.sel.contours[[popmats[i,2]]],cont=c(50,80,90,95),col=contcols2[5:2],add=T,drawpoints=F,pch=16,cex=.5,col.pt=rgb(0,1,0,alpha=.1))
	plot(all.sel.contours[[popmats[i,3]]],cont=c(50,80,90,95),col=contcols3[5:2],add=T,drawpoints=F,pch=16,cex=.5,col.pt=rgb(1,0,1,alpha=.1))
	mtext(popnamestemp[i],side=3, line=1, cex=1.2)
	abline(h=0,lty=2);abline(v=0,lty=2)
#	if(i==1){mtext("Selection for later flowering",side=3,line=2,cex=1.4)#}
	if(i==5){mtext("Root Mass Response",side=2,line=3,cex=1.5)}
	if(i==9){mtext("Flowering Time Response",side=1,line=4,cex=1.5,adj=-8)}
	legend(-3.75,.75,c("Biota15.0","Biota14.3","Biota13.0"),col=c(rgb(0,0,0),rgb(0,1,0),rgb(0,0,1)),pch=16,bty="n")
}
dev.off()
