#makes heritability figure
library(MCMCglmm)

load("~/GenoPheno/ComparisonG/popsoil.animal.mods.5trait.Rdata") #G matrices object: popsoil5.animal.mods


####G matrices
AnimalModelheritability <-function(mod,gcols,ecols,ntraits){
hmcmc <- mod$VCV[,gcols]/ (mod$VCV[,gcols]+mod$VCV[,ecols])
hmean <- diag(matrix(colMeans(hmcmc),ncol=ntraits))#genetic (co)variance for traits
SEh <- array(dim=c(ntraits,2))
SEh[,1] <- diag(matrix(HPDinterval(hmcmc)[,1],ncol=ntraits))
SEh[,2] <- diag(matrix(HPDinterval(hmcmc)[,2],ncol=ntraits))#upper and lower 95% interval of posterior for heritability
return(list(hsqr= hmean,interval= SEh))
}

allpopsoilProcess <- lapply(popsoil5.animal.mods, function(z) AnimalModelheritability(z,gcols=1:25,ecols=26:50,ntraits=5))

Gnames.ps <- paste(rep(c("cl","cu","da","fp","m","mc","mt","tc","tx","tz"),times=3),".",rep(c("mt","mc","tc"),each=10),sep="")#represents actual order of Gs. 
#Gnames.ps useful for checking rearranging vectors to make sure they rearrange things as you think.
#rearranging vector (mean annual temp sort):
swappopsinsoil <- rep(c(10,5,9,7,3,4,6,1,8,2),each=3) + rep(c(0,10,20),times=10)#note funny order of 0,10,20 retaining temp order of soils which was happenstance of G model list

hpop_ordered <- lapply(1:30, function(z) allpopsoilProcess[[swappopsinsoil[z]]])

unlist_h <- lapply(1:5, function(k)  unlist(lapply(hpop_ordered, function(z) z[[1]][k]))) #all means
unlist_lo <- lapply(1:5, function(k)  unlist(lapply(hpop_ordered, function(z) z[[2]][k,1]))) #all lower interval bound
unlist_up <- lapply(1:5, function(k)  unlist(lapply(hpop_ordered, function(z) z[[2]][k,2]))) #upper bound

traitnames <- c("Days to flowering","Days to germination","Shoot mass","Root Mass","Stem width")
herit <- c("","","heritability","","")

pdf("~/GenoPheno/ComparisonG/heritabilities.pdf",height=6,width=4)
par(mar=c(2,4,2,1))
par(mfrow=c(5,1))
par(oma=c(2,0,0,0))
for(i in 1:5){
plot(unlist_h[[i]]~c(1:30),ylim=c(0,1),pch=16,col=c(rgb(0.8,0,0),rgb(.5,0,.5),rgb(0,0,0.8)),xaxt="n",
	xlab="",ylab="",yaxt="n",main="")
mtext(traitnames[i],side=3,line=0.2)
mtext(herit[i],side=2,line=2.2)
axis(side=1,at=c(2,5,8,11,14,17,20,23,26,29),labels=c("19.8","18.6","15.3","15.0","14.7","14.4","14.3","13.2","13.0","12.9"))#c("TZ","ML","TX","MT","DA","FP","MC","CL","TC","CU")
axis(side=2,at=c(0,1),labels=c("0","1"))
abline(v=seq(from=3.5, to=27.5,by=3),lty=3,col=rgb(0,0,0,.5))
abline(h=0.1,lty=2)
arrows(1:30,unlist_lo[[i]],1:30,unlist_up[[i]],length=0, col = c(rgb(0.8,0,0),rgb(.5,0,.5),rgb(0,0,0.8)))
}
mtext("Plant population source MAT",side=1,line=2.5)
dev.off()