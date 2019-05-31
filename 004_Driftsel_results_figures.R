library(driftsel)
library(SDMTools) #legend.gradient

load('~/Dropbox/AO-3 mine/scaled traits/scaledtraits.samp.mc_med_LG95_440k.Rdata')
load('~/Dropbox/AO-3 mine/scaled traits/scaledtraits.samp.mt_med_LG95_440k.Rdata')
load('~/Dropbox/AO-3 mine/scaled traits/scaledtraits.samp.tc_med_LG95_440k.Rdata')

load('~/Dropbox/AO-3 mine/scaled traits/scaledtraits.samp.mc.h_med_LG95_440k.Rdata')
load('~/Dropbox/AO-3 mine/scaled traits/scaledtraits.samp.mt.h_med_LG95_440k.Rdata')
load('~/Dropbox/AO-3 mine/scaled traits/scaledtraits.samp.tc.h_med_LG95_440k.Rdata')

#useful in color plotting
range01=function(x){
newnums=(x-min(x,na.rm=T))/(max(x,na.rm=T)-min(x,na.rm=T))
return(newnums)
}

##MODIFIED FUNCTIONS FROM PACKAGE DRIFTSEL, code is nearly identical to the driftsel function without the "mod." in the name
mod.ellipse <- function (mu, Sig, prob.mass, lwd=1, col="black") 
{
    radius <- sqrt(qchisq(prob.mass, 2))
    theta <- 2 * pi * (0:360)/360
    unit.circle <- cbind(cos(theta), sin(theta))
    Q <- base::chol(Sig, pivot = TRUE)
    order <- order(attr(Q, "pivot"))
    ellipse <- t(mu + radius * t(unit.circle %*% Q[, order]))
    lines(ellipse[, 1], ellipse[, 2], lwd = lwd, col = col,lty=3)
}
mod.viz.theta <- function (ds.out, traits,colvec, siz = 0.5, main = NA, optylim = NULL,optxlim=NULL,ylab=NULL,xlab=NULL, pltpts = TRUE,cex=1,cex.axis=1,cex.main=1) {
	popef = apply(ds.out$pop.ef[, traits, ],c(1,2),mean)
    G = apply(ds.out$G[traits, traits, ],c(1,2),mean)
    TH = apply(ds.out$theta,c(1,2),mean)
    mu = apply(ds.out$fixed.ef[1, traits, ],c(1),mean)#,2? in c(1)
    npop = ncol(TH)
    i = 1
    j = 2
    ei = mu[i]
    ej = mu[j]
    Gthis = matrix(c(G[i, i], G[i, j], G[j, i], G[j, j]), ncol = 2)
    xlab = ifelse(is.character(xlab),xlab, paste("trait", traits[i]))
    ylab = ifelse(is.character(ylab),ylab, paste("trait", traits[j]))
    M = max(c(sqrt(2 * G[i, i] * diag(TH)), sqrt(2 * G[j, j] * 
        diag(TH)), popef))
    if(is.vector(optylim) & is.vector(optxlim)){plot(ei, ej, pch = 16, cex = 0, xlim = optxlim, 
        ylim = optylim, xlab = xlab, ylab = ylab, main = main,cex.axis=cex.axis,cex.main=cex.main)
    } else {plot(ei, ej, pch = 16, cex = 0, xlim =  c(ei - M, ei + M), 
        ylim = c(ej - M, ej + M), xlab = xlab, ylab = ylab, main = main)}
    for (k in 1:npop) {
        mu = c(ei, ej)
        Sigma = 2 * TH[k, k] * Gthis
   #     if (k > 9) {
   #         k = 1 + k%%8
   #     }
        mod.ellipse(mu, Sigma, siz, 1, colvec[k])
    }
    if(pltpts==TRUE){
    points(ei + popef[, i], ej + popef[, j], cex = cex, col = colvec,pch=16)}
    points(ei, ej, pch=16,cex=1.5)
}


#which populations and traits fall outside one-dimensional intervals?
getoutside <- function(dsobj,ntrait, npop,CI=.95){
tunit <- cbind(cos(2 * pi * (0:360)/360), sin(2 * pi * (0:360)/360))
trad <- sqrt(qchisq(CI, 2))
popef <- apply(dsobj$pop.ef,c(1,2),mean)
fixefs <- apply(dsobj$fixed.ef,c(1,2),mean)
popoutside <- matrix(NA, ncol=ntrait,nrow=npop)
popmeandiff <- matrix(NA, ncol=ntrait,nrow=npop)
th <- apply(dsobj$theta,c(1,2),mean)
for(k in 1:npop){
	ranges <- matrix(NA,ncol=2,nrow=ntrait)
	for(i in 1:(ntrait-1)){
		tg <- apply(dsobj$G,c(1,2),mean)[c(i,i+1),c(i,i+1)]
		tfx <- apply(dsobj$fixed.ef,c(1,2),mean)[c(i,i+1)]
		tq <-  base::chol(2 * th[i, i] * tg, pivot = TRUE)
		oq <- order(attr(tq, "pivot"))
		tempE <- t(tfx + trad * t(tunit %*% tq[, oq]))#this is the ellipse from mod.ellipse and mod.viztheta
		if(i < 2){ ranges[i,]<- range(tempE[,1])}
		ranges[i+1,]<- range(tempE[,2])#here we just project it to the trait axes
		}
	popoutside[k,] <-  (popef[k,]+fixefs) > ranges[,2] | (popef[k,]+fixefs) < ranges[,1]
	popmeandiff[k,] <- popef[k,]
	}
return(list(popoutside, popmeandiff))
}

# using a polygon approximation of the ellipse, i.e. which populations and traits fall outside expectations in two-dimensional traits
getoutsideBiVar <- function(dsobj,ntrait, npop,CI=.95){
tunit <- cbind(cos(2 * pi * (0:360)/360), sin(2 * pi * (0:360)/360))
trad <- sqrt(qchisq(CI, 2))
popef <- apply(dsobj$pop.ef,c(1,2),mean)
fixefs <- apply(dsobj$fixed.ef,c(1,2),mean)
popoutside <- matrix(NA, ncol=ntrait,nrow=npop)
#popmeandiff <- matrix(NA, ncol=ntrait,nrow=npop)
th <- apply(dsobj$theta,c(1,2),mean)
outmats <- list()
for(k in 1:npop){
	isout <- matrix(NA,ncol=ntrait,nrow=ntrait) #need to fill a number of cells equiv to (1:(ntrait-1))
	for(i in 1:(ntrait-1)){
		for(j in (i+1):(ntrait)){
		tg <- apply(dsobj$G,c(1,2),mean)[c(i,j),c(i,j)]
		tfx <- apply(dsobj$fixed.ef,c(1,2),mean)[c(i,j)]
		tq <-  base::chol(2 * th[i, i] * tg, pivot = TRUE)
		oq <- order(attr(tq, "pivot"))
		tempE <- t(tfx + trad * t(tunit %*% tq[, oq]))
		tpp <- cbind(popef[k,i]+fixefs[i],popef[k,j]+fixefs[j]) #am assuming the x y relationship here in the ellipse points, baseed on what went in.
		dists <- sqrt( (tpp[,1] - tempE[,1])^2 + (tpp[,2] - tempE[,2])^2)
		mincoord <- tempE[which(dists == min(dists))[1],]
		distCenterfromE <- sqrt((mincoord[1] - fixefs[i])^2 + (mincoord[2] - fixefs[j])^2)
		distCenterfromPoint <- sqrt((tpp[1] - fixefs[i])^2 + (tpp[2] - fixefs[j])^2)
		isout[i,j] <- distCenterfromPoint > distCenterfromE
		}
		
		}
	outmats[[k]] <- isout
	}
	popout <- matrix(0,ncol=ntrait,nrow=ntrait)
	for(k in 1:npop){popout <- popout + outmats[[k]]}	
return(list(outmats,popout))
}



#################################################

#############FIGURES

envmatall <-read.csv("~/Dropbox/AO-3 mine/AlphabeticalPopEnvDat.csv")
envmat <- envmatall[,c(1,4,48)]#just mean temp, annual precip, and soil water holding capacity

# color palettes
wtob <- colorRampPalette(c(rgb(1,1,1),rgb(0,0,0.5)))
rb <-colorRampPalette(c(rgb(0,0,1),rgb(1,0,0)))
gp <-colorRampPalette(c(rgb(0.5,0,0.75),rgb(1,1,1),rgb(0,0.5,0)))


alphtomat <- order(envmat$TAnn,decreasing=T)

tempnames <- c("19.8","18.6","15.3","15.0","14.7","14.4","14.3","13.2","13.0","12.9")
prettytraitnames <- c("Days to flowering","Days to germination","Tassel length","Shoot biomass","Root biomass","Height","Stem width","Leaf width","Stem greenness")


pdf("~/Dropbox/AO-3 mine/scaled traits/dtf and stw 440k.pdf",width=8,height=3)
layout(matrix(1:4,ncol=4),widths=c(3,3,3,2.25))
#par(mfrow=c(1,3))
par(mar=c(6,5,5,1))
popcolvec  <- 	rgb(range01(envmat$TAnn),0,1-range01(envmat$TAnn))
mod.viz.theta(samp.mt,c(1,7),siz=.95,colvec=popcolvec,optylim=c(-1,2),optxlim=c(-1,2),xlab="",ylab="",main="Biota15.0",cex=1,cex.axis=1.5,cex.main=2)
mtext("Scaled Stem Width",side=2, line=2.5,cex=1.2)
mod.viz.theta(samp.mc,c(1,7),siz=.95,colvec=popcolvec,optylim=c(-1,2),optxlim=c(-1,2),ylab="",xlab="",main="Biota14.3",cex=1,cex.axis=1.5,cex.main=2)
mtext("Scaled Days to Flowering",side=1, line=3,cex=1.2)
mod.viz.theta(samp.tc,c(1,7),siz=.95,colvec=popcolvec,optylim=c(-1,2),optxlim=c(-1,2),xlab="",ylab="",main="Biota13.0",cex=1,cex.axis=1.5,cex.main=2)
par(mar=c(5,6,3,5))
image(1,seq(from=.5,to=9.5,by=1),t(matrix(rev(sort(envmat$TAnn)),nrow=10)),col=rb(50),zlim=c(129,198),ylim=c(0,10),ylab="",xlab="",bty="n",xaxt="n",yaxt="n")
#par(mfrow=c(1,3))
axis(side=2, at=seq(from=.5,to=9.5,by=1), labels=tempnames, las=2, cex.axis = 1.5)
mtext("Plant source MAT", side=3, line =1, cex =1)
dev.off()


poptr.mt <- getoutside(samp.mt,9,10,CI=.95)
poptr.mc <- getoutside(samp.mc,9,10,CI=.95)
poptr.tc <- getoutside(samp.tc,9,10,CI=.95)


poptr.mt.9 <- getoutside(samp.mt,9,10,CI=.9)
poptr.mc.9 <- getoutside(samp.mc,9,10,CI=.9)
poptr.tc.9 <- getoutside(samp.tc,9,10,CI=.9)




divrange <-  round(c(-1*max( max(abs(poptr.mt[[2]])),max(abs(poptr.mc[[2]])),max(abs(poptr.tc[[2]])) ),max( max(abs(poptr.mt[[2]])),max(abs(poptr.mc[[2]])),max(abs(poptr.tc[[2]])) ) ),digits=1)
#based on this, 
pdf("~/Dropbox/AO-3 mine/scaled traits/sigdiffs by pop and trait CI 95 names means torder 440.pdf",width=8,height=4)
	layout(matrix(1:5,ncol=5),widths=c(1.6,3,3,3,1.5))
	par(mar=c(13,6,3,0))
	image(1,seq(from=.5,to=9.5,by=1),t(matrix(rev(sort(envmat$TAnn)),nrow=10)),col=rb(50),zlim=c(129,198),ylim=c(0,10),ylab="",xlab="",bty="n",xaxt="n",yaxt="n")
		axis(side=2, at=seq(from=.5,to=9.5,by=1), labels=tempnames, las=2, cex.axis = 1.5)
		mtext("Plant population MAT (row)", side=2, line =3.75, cex =1)
		par(mar=c(13,1,3,1))
	image(1:9,seq(from=.5,to=9.5,by=1),t(poptr.mt[[2]][alphtomat,]),ylim=c(0,10),axes=F,ylab="",xlab="",col=gp(50),zlim=divrange,main="Biota15.0",cex.main=1.5)#Mateo tezoquipan, mat 15.0
		points(0,0)
		for(i in 1:10){
			yval <- seq(from=.5,to=9.5,by=1)[i]
			points(rep(yval,times=9)~c(1:9), pch= ifelse(as.numeric(poptr.mt[[1]][alphtomat,][i,])>0, "*",NA),cex=2)
			points(rep(yval,times=9)~c(1:9), pch= ifelse(as.numeric(poptr.mt.9[[1]][alphtomat,][i,]) - as.numeric(poptr.mt[[1]][alphtomat,][i,]) >0, 1,NA),cex=1)
		}
		axis(side=1, at=1:9, labels=prettytraitnames, las=2, cex.axis = 1.5)
	par(mar=c(13,1,3,1))
	image(1:9,seq(from=.5,to=9.5,by=1),t(poptr.mc[[2]][alphtomat,]),ylim=c(0,10),axes=F,ylab="",xlab="",col=gp(50),zlim=divrange,main="Biota14.3",cex.main=1.5)#matias cuijingo, mat 14.3
		points(0,0)
		for(i in 1:10){
			yval <- seq(from=.5,to=9.5,by=1)[i]
			points(rep(yval,times=9)~c(1:9), pch= ifelse(as.numeric(poptr.mc[[1]][alphtomat,][i,])>0, "*",NA),cex=2)
			points(rep(yval,times=9)~c(1:9), pch= ifelse(as.numeric(poptr.mc.9[[1]][alphtomat,][i,]) - as.numeric(poptr.mc[[1]][alphtomat,][i,])>0, 1,NA),cex=1)
		}
		axis(side=1, at=1:9, labels=prettytraitnames, las=2, cex.axis = 1.5)
	par(mar=c(13,1,3,1))
	image(1:9,seq(from=.5,to=9.5,by=1),t(poptr.tc[[2]][alphtomat,]),ylim=c(0,10),axes=F,ylab="",xlab="",col=gp(50),zlim=divrange,main="Biota13.0",cex.main=1.5)#toluca C, mat 14.3
		points(0,0)
		for(i in 1:10){
			yval <- seq(from=.5,to=9.5,by=1)[i]
			points(rep(yval,times=9)~c(1:9), pch= ifelse(as.numeric(poptr.tc[[1]][alphtomat,][i,])>0, "*",NA),cex=2)
			points(rep(yval,times=9)~c(1:9), pch= ifelse(as.numeric(poptr.tc.9[[1]][alphtomat,][i,])-as.numeric(poptr.tc[[1]][alphtomat,][i,])>0, 1,NA),cex=1)
		}
		axis(side=1, at=1:9, labels=prettytraitnames, las=2, cex.axis = 1.5)
	par(mar=c(13,0,1,0))
	plot(c(1,10)~c(1,2),xaxt="n", yaxt="n",bty="n",pch=NA,ylab="",xlab="")
		legend.gradient(cbind(c(1,1.25,1.25,1),c(1,1,9,9)),cols=gp(50),title="divergence",limits=as.character(divrange),cex=1.5)
dev.off()


bivar.mt <- getoutsideBiVar(samp.mt,9,10,CI=.95)
bivar.mc <- getoutsideBiVar(samp.mc,9,10,CI=.95)
bivar.tc <- getoutsideBiVar(samp.tc,9,10,CI=.95)
WtoB <- colorRampPalette(c(rgb(.9,.9,.9),rgb(0,0,1),rgb(0,0,.25)))
WtoW <- colorRampPalette(c(rgb(1,1,1),rgb(1,1,1)))

#number of traits in any trait combo
# note bc of flowering time univariate always least one pop always detected so any trait in combo with flowering time must have at least one pop outside
#number of pops detected div in any trait combo, (out of 10 poss per soil)
sum(unlist( lapply(1:length(bivar.mt[[1]]),function(z) sign(sum(unlist(bivar.mt[[1]][z]),na.rm=T )) )))
sum(unlist( lapply(1:length(bivar.mc[[1]]),function(z) sign(sum(unlist(bivar.mc[[1]][z]),na.rm=T )) )))
sum(unlist( lapply(1:length(bivar.tc[[1]]),function(z) sign(sum(unlist(bivar.tc[[1]][z]),na.rm=T )) )))
#total detections (out of 360 poss per soil)
sum(unlist( lapply(1:length(bivar.mt[[1]]),function(z) (sum(unlist(bivar.mt[[1]][z]),na.rm=T )) ) ))
sum(unlist( lapply(1:length(bivar.mc[[1]]),function(z) (sum(unlist(bivar.mc[[1]][z]),na.rm=T )) ) ))
sum(unlist( lapply(1:length(bivar.tc[[1]]),function(z) (sum(unlist(bivar.tc[[1]][z]),na.rm=T )) ) ))
#times each population detected outside as a vector
bivar.mt.pops <- unlist( lapply(1:length(bivar.mt[[1]]),function(z) sum(unlist(bivar.mt[[1]][z]),na.rm=T ) ) ) 
bivar.mc.pops <- unlist( lapply(1:length(bivar.mc[[1]]),function(z) sum(unlist(bivar.mc[[1]][z]),na.rm=T ) ) )
bivar.tc.pops <- unlist( lapply(1:length(bivar.tc[[1]]),function(z) sum(unlist(bivar.tc[[1]][z]),na.rm=T ) ) )


wtog2 <- colorRampPalette(c(rgb(0.8,0.8,0.8),rgb(0,1,0)) )

pdf("~/Dropbox/AO-3 mine/scaled traits/bivariate sigdiffs by trait CI 95 names means torder 440.pdf",width=13,height=7)
 layout(matrix(c(1:21),nrow=3,ncol=7),widths=c(rep(c(0.4,2),times=3),0.3),heights=rev(c(0.5,1.25,2)))
 par(oma = c(5,12,2,3))
 #mt
  par(mar=c(1,1,1,1))
 image( matrix((rev(colSums(poptr.mt[[1]]))),nrow=1),col=WtoB(11),xaxt="n",yaxt="n",zlim=c(0,10))
axis(side=2, at=seq(from=0,to=1,length.out=9), labels=rev(prettytraitnames), las=1, cex.axis = 1.5)
 image( matrix((colSums(poptr.mt[[1]])),nrow=1),col=WtoW(3),xaxt="n",yaxt="n",bty="n") 
 image( matrix((colSums(poptr.mt[[1]])),nrow=1),col=WtoW(3),xaxt="n",yaxt="n",bty="n") 
 image( bivar.mt[[2]][,9:1],col=WtoB(11),xaxt="n",yaxt="n",zlim=c(0,10),main="Biota15.0",cex.main=1.5)
  par(mar=c(12,1,1,1))
 image( matrix((colSums(poptr.mt[[1]])),ncol=1),col=WtoB(11),xaxt="n",yaxt="n",zlim=c(0,10))
 axis(side=1, at=seq(from=0,to=1,length.out=9), labels=prettytraitnames, las=3, cex.axis = 1.5)
  par(mar=c(1,1,1,1))
image( matrix(bivar.mt.pops[alphtomat],ncol=1),col=wtog2(35),xaxt="n",yaxt="n",zlim=c(1,36))
 axis(side=1, at=seq(from=0,to=1,length.out=10), labels=tempnames, las=3, cex.axis = 1.5)
axis(side=2, at=c(0), labels=c("Contributing populations"), las=1, cex.axis = 1.5)
	mtext("MAT of population",side = 1,line=4)
#mc
 image( matrix((rev(colSums(poptr.mc[[1]]))),nrow=1),col=WtoB(11),xaxt="n",yaxt="n",zlim=c(0,10))
 image( matrix((colSums(poptr.mc[[1]])),nrow=1),col=WtoW(3),xaxt="n",yaxt="n",zlim=c(0,10),bty="n")
 image( matrix((colSums(poptr.mc[[1]])),nrow=1),col=WtoW(3),xaxt="n",yaxt="n",zlim=c(0,10),bty="n")
 image( bivar.mc[[2]][,9:1],col=WtoB(11),xaxt="n",yaxt="n",main="Biota14.3",cex.main=1.5,zlim=c(0,10))
   par(mar=c(12,1,1,1))
image( matrix((colSums(poptr.mc[[1]])),ncol=1),col=WtoB(11),xaxt="n",yaxt="n",zlim=c(0,10))
 axis(side=1, at=seq(from=0,to=1,length.out=9), labels=prettytraitnames, las=3, cex.axis = 1.5)
  par(mar=c(1,1,1,1))
 image( matrix(bivar.mc.pops[alphtomat],ncol=1),col=wtog2(35),xaxt="n",yaxt="n",zlim=c(1,36))
 axis(side=1, at=seq(from=0,to=1,length.out=10), labels=tempnames, las=3, cex.axis = 1.5)
	mtext("MAT of population",side = 1,line=4)
#tc
 image( matrix((rev(colSums(poptr.tc[[1]]))),nrow=1),col=WtoB(11),xaxt="n",yaxt="n",zlim=c(0,10))
 image( matrix((colSums(poptr.tc[[1]])),nrow=1),col=WtoW(3),xaxt="n",yaxt="n",zlim=c(0,10),bty="n")
 image( matrix((colSums(poptr.tc[[1]])),nrow=1),col=WtoW(3),xaxt="n",yaxt="n",zlim=c(0,10),bty="n")
 image( bivar.tc[[2]][,9:1],col=WtoB(11),xaxt="n",yaxt="n",main="Biota13.0",cex.main=1.5,zlim=c(0,10))
   par(mar=c(12,1,1,1))
image( matrix((colSums(poptr.tc[[1]])),ncol=1),col=WtoB(11),xaxt="n",yaxt="n",zlim=c(0,10))
  axis(side=1, at=seq(from=0,to=1,length.out=9), labels=prettytraitnames, las=3, cex.axis = 1.5)
   par(mar=c(1,1,1,1))
 image( matrix(bivar.tc.pops[alphtomat],ncol=1),col=wtog2(35),xaxt="n",yaxt="n",zlim=c(1,36))
 axis(side=1, at=seq(from=0,to=1,length.out=10), labels=tempnames, las=3, cex.axis = 1.5)
	mtext("MAT of population",side = 1,line=4)
#legend
  par(mar=c(1,2,1,1))
 image( matrix(0:10,nrow=1),col=WtoB(11),xaxt="n",yaxt="n",zlim=c(0,10))
axis(side=4, at=c(0,1), labels=c("0","10"), las=1, cex.axis = 1.5)
mtext("# significant populations",side=2)
 image( matrix(1:36,nrow=1),col=wtog2(35),xaxt="n",yaxt="n",zlim=c(1,36))
axis(side=4, at=c(0,1), labels=c("1","36"), las=1, cex.axis = 1.5)
mtext("# significant bivariate cases",side=2)
 image( matrix(0:36,nrow=1),col=WtoW(36),xaxt="n",yaxt="n",zlim=c(1,36),bty="n")
dev.off()
