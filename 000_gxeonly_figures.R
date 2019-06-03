
range01=function(x){
newnums <- sapply(1:length(x), function(z) ifelse(!is.na(x[z]), (x[z]-min(x,na.rm=T))/(max(x,na.rm=T)-min(x,na.rm=T)), NA ))
return(newnums)
}


####READ IN DATA
soils <- read.csv("~/GenoPheno/covmat.csv",stringsAsFactors=F,header=T)
allpheno <- read.csv("~/GenoPheno/phenomat.csv",stringsAsFactors=F,header=T)
allped <- read.csv("~/GenoPheno/pedmat.csv",stringsAsFactors=F,header=T)
Apopenvdat <-read.csv("AlphabeticalPopEnvDat.csv",stringsAsFactors=F,header=T)

greens <- (allpheno$green)/(allpheno$red + allpheno$blue + allpheno$green)
#swaps colors to ratios with total brightness
allpheno <- allpheno[,-c(3,10:12,14,16)]#cuts some traits not analyzed for various reasons (some colinear, some non-normal, some invariant), leaves 9 traits
allpheno$green <- greens
allpheno$germday <- allpheno$germday - 16 # is in days since july 1st, and later want it in Days to Germination, so here the script subtracts planting date

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

#### FIGURE
##ploting pop-momxe effects.
#plot phenotypes across soil treatments
# color by pop source temp. plot all individuals connected by lines was original. new figure plots differently.
traits.mt <-allpheno[which(soils$treatment==3),]
traits.mc <- allpheno[which(soils$treatment==5),]
traits.tc <- allpheno[which(soils$treatment==7),]

ped.mt <- allped[which(soils$treatment==3),] 
ped.mc <- allped[which(soils$treatment==5),]
ped.tc <- allped[which(soils$treatment==7),]

damtrait.mc <- sapply(2:10,function(z) tapply(traits.mc[,z],ped.mc$dam,mean,na.rm=T))
damtrait.mt <- sapply(2:10,function(z) tapply(traits.mt[,z],ped.mt$dam,mean,na.rm=T))
damtrait.tc <- sapply(2:10,function(z) tapply(traits.tc[,z],ped.tc$dam,mean,na.rm=T))
damtraitse.mc <- sapply(2:10,function(z) tapply(traits.mc[,z],ped.mc$dam,sd,na.rm=T)) / sqrt(sapply(2:10,function(z) tapply(sign(traits.mc[,z]),ped.mc$dam,sum,na.rm=T)))
damtraitse.mt <- sapply(2:10,function(z) tapply(traits.mt[,z],ped.mt$dam,sd,na.rm=T)) / sqrt(sapply(2:10,function(z) tapply(sign(traits.mt[,z]),ped.mt$dam,sum,na.rm=T)))
damtraitse.tc <- sapply(2:10,function(z) tapply(traits.tc[,z],ped.tc$dam,sd,na.rm=T)) / sqrt(sapply(2:10,function(z) tapply(sign(traits.tc[,z]),ped.tc$dam,sum,na.rm=T)))


poptrait.mc <- sapply(2:10,function(z) tapply(traits.mc[,z],ped.mc$dam.pop,mean,na.rm=T))
poptrait.mt <- sapply(2:10,function(z) tapply(traits.mt[,z],ped.mt$dam.pop,mean,na.rm=T))
poptrait.tc <- sapply(2:10,function(z) tapply(traits.tc[,z],ped.tc$dam.pop,mean,na.rm=T))
poptraitse.mc <- sapply(2:10,function(z) tapply(traits.mc[,z],ped.mc$dam.pop,sd,na.rm=T)) / sqrt(sapply(2:10,function(z) table(ped.mc$dam.pop[which( !is.na(traits.mc[,z]) )])))
poptraitse.mt <- sapply(2:10,function(z) tapply(traits.mt[,z],ped.mt$dam.pop,sd,na.rm=T)) / sqrt(sapply(2:10,function(z) table(ped.mc$dam.pop[which( !is.na(traits.mt[,z]) )])))
poptraitse.tc <- sapply(2:10,function(z) tapply(traits.tc[,z],ped.tc$dam.pop,sd,na.rm=T)) / sqrt(sapply(2:10,function(z) table(ped.mc$dam.pop[which( !is.na(traits.tc[,z]) )])))


damtraits.soils <- array(,dim=c(100,3,9))#100 moms, 3 soils, 9 traits ##REORDERED SOILS
for(i in 1:9){
	damtraits.soils[,2,i] <- damtrait.mc[,i]
	damtraits.soils[,1,i] <- damtrait.mt[,i]
	damtraits.soils[,3,i] <- damtrait.tc[,i]
}
damtraitsse.soils <- array(,dim=c(100,3,9))#100 moms, 3 soils, 9 traits ##REORDERED SOILS
for(i in 1:9){
	damtraitsse.soils[,2,i] <- damtraitse.mc[,i]
	damtraitsse.soils[,1,i] <- damtraitse.mt[,i]
	damtraitsse.soils[,3,i] <- damtraitse.tc[,i]
}

poptraits.soils <- array(,dim=c(10,3,9))#10 pops, 3 soils, 9 traits ## REORDERED SOILS
for(i in 1:9){
	poptraits.soils[,2,i] <- poptrait.mc[,i]
	poptraits.soils[,1,i] <- poptrait.mt[,i]
	poptraits.soils[,3,i] <- poptrait.tc[,i]
}
poptraitsse.soils <- array(,dim=c(10,3,9))#10 pops, 3 soils, 9 traits ## REORDERED SOILS
for(i in 1:9){
	poptraitsse.soils[,2,i] <- poptraitse.mc[,i]
	poptraitsse.soils[,1,i] <- poptraitse.mt[,i]
	poptraitsse.soils[,3,i] <- poptraitse.tc[,i]
}

tc.mn <- apply(traits.tc[,2:10],2,mean,na.rm=T)
tc.se <- apply(traits.tc[,2:10],2,sd,na.rm=T)/sqrt(colSums(!is.na(traits.tc[,2:10])))
mc.mn <- apply(traits.mc[,2:10],2,mean,na.rm=T)
mc.se <- apply(traits.mc[,2:10],2,sd,na.rm=T)/sqrt(colSums(!is.na(traits.mc[,2:10])))
mt.mn <- apply(traits.mt[,2:10],2,mean,na.rm=T)
mt.se <- apply(traits.mt[,2:10],2,sd,na.rm=T)/sqrt(colSums(!is.na(traits.mt[,2:10])))

traitcols <- c("Days to flowering","Days to germination","Tassel length cm","Shoot biomass g","Root biomass g","Height cm","Stem width mm","Leaf width mm","Stem greenness %")#nice names
rbb<-rgb(range01(Apopenvdat$TAnn),0,1-range01(Apopenvdat$TAnn))
rbsoft <- rgb(range01(Apopenvdat$TAnn),0,1-range01(Apopenvdat$TAnn),alpha=.75)
rbsoft1<- rgb(range01(Apopenvdat$TAnn),0,1-range01(Apopenvdat$TAnn),alpha=.3)
rbsoft2<- rgb(range01(Apopenvdat$TAnn),0,1-range01(Apopenvdat$TAnn),alpha=.1)
damtoalphpop <- tapply(allped$dam.pop,allped$dam,mean)


pdf("~/popdamtraits_mompopsepSE_soilpan_3onlycolY-axis.pdf",height=6,width=3)
	layout(matrix(1:9,ncol=3,byrow=T),heights=c(5,5,5),widths=c(0.75,0.9,1.1))
		jitter1 <- jitter(rep(1,times=10),factor=3)
		jitter2 <- jitter(rep(2,times=10),factor=3)
		jitter3 <- jitter(rep(3,times=10),factor=3)
		jitter1d <- jitter(rep(1,times=100),factor=3)
		jitter2d <- jitter(rep(2,times=100),factor=3)
		jitter3d <- jitter(rep(3,times=100),factor=3)
	par(oma = c(0,2,2,0))
	for(i in c(5,3,1)){
	par(mar=c(2,3,1,0))
		ylims <- c(min(c(mt.mn[i]-mt.se[i],mc.mn[i]-mc.se[i],tc.mn[i]-tc.se[i])),max(c(mt.mn[i]+mt.se[i],mc.mn[i]+mc.se[i],tc.mn[i]+tc.se[i])) )
		ylimp <- c(min(poptraits.soils[,,i]-poptraitsse.soils[,,i]),max(poptraits.soils[,,i]+poptraitsse.soils[,,i]))
		plot(c(mt.mn[i],mc.mn[i],tc.mn[i])~c(1,2,3),pch=1,cex=1,ylab="",xlab="",xaxt="n",cex.lab=1,xlim=c(0.5,3.5), bty="n",
				ylim= ylimp)
		axis(seq(from=-1000,to=1000,length.out=2),side=2, labels=NULL)
		arrows(c(1,2,3),c(mt.mn[i]-mt.se[i],mc.mn[i]-mc.se[i],tc.mn[i]-tc.se[i]),y1=c(mt.mn[i]+mt.se[i],mc.mn[i]+mc.se[i],tc.mn[i]+tc.se[i]),length=0  )#,add=T)
		mtext(traitcols[i],side=2, line=2)
	par(mar=c(2,0,1,1))
		plot(poptraits.soils[,,i]~1, pch=NA,ylab="",xlab="",xaxt="n",yaxt="n",cex.lab=1,xlim=c(0.8,3.2), bty = "n", ylim=ylimp)
		rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = rgb(0,0,0,alpha=0.1),lty=NULL,border=NA)
		sapply(1:nrow(poptraits.soils[,,i]), function(z) lines(c(jitter1[z],jitter2[z],jitter3[z]),poptraits.soils[z,,i],col=rbsoft[z],lty=1 ))
		arrows(jitter1,poptrait.mt[,i]-poptraitse.mt[,i],y1=poptrait.mt[,i]+poptraitse.mt[,i],col=rbsoft,length=0  )#,add=T)
		arrows(jitter2,poptrait.mc[,i]-poptraitse.mc[,i],y1=poptrait.mc[,i]+poptraitse.mc[,i],col=rbsoft,length=0  )#,add=T)
		arrows(jitter3,poptrait.tc[,i]-poptraitse.tc[,i],y1=poptrait.tc[,i]+poptraitse.tc[,i],col=rbsoft,length=0  )#,add=T)
	par(mar=c(2,1,1,1))
		ylimf <- c( min( c(min(damtraits.soils[,,i]-damtraitsse.soils[,,i],na.rm=T),min(damtraits.soils[,,i],na.rm=T)) ),
			max( c(max(damtraits.soils[,,i]+damtraitsse.soils[,,i],na.rm=T),max(damtraits.soils[,,i],na.rm=T)) ) )
		plot(damtraits.soils[,,i]~1, pch=NA,ylab="",xlab="",xaxt="n",cex.lab=1,xlim=c(0.8,3.2), bty="n", ylim=ylimf )
		axis(seq(from=-1000,to=1000,length.out=2),side=2, labels=NULL)
		sapply(1:nrow(damtraits.soils[,,i]), function(z) lines(c(jitter1d[z],jitter2d[z],jitter3d[z]),damtraits.soils[z,,i],col=rbsoft1[damtoalphpop[z]],lty=1 ))
		arrows(jitter1d,damtrait.mt[,i]-damtraitse.mt[,i],y1=damtrait.mt[,i]+damtraitse.mt[,i],col=rbsoft1[damtoalphpop],length=0  )#,add=T)
		arrows(jitter2d,damtrait.mc[,i]-damtraitse.mc[,i],y1=damtrait.mc[,i]+damtraitse.mc[,i],col=rbsoft1[damtoalphpop],length=0  )#,add=T)
		arrows(jitter3d,damtrait.tc[,i]-damtraitse.tc[,i],y1=damtrait.tc[,i]+damtraitse.tc[,i],col=rbsoft1[damtoalphpop],length=0  )#,add=T)
	}
dev.off()


# 
pdf("~/popdamtraits_soilsOct2018_optionsMomPOPSOIL.pdf",height=8,width=8)
#all traits by mom,  soils on different lines
	layout(matrix(1:27,ncol=3,byrow=T),width=c(4,1,0.75))
	par(oma=c(0,0,1,2))
	par(mar=c(2,3,0,0))
	for(i in 1:9){
		trtmns <- as.vector(t(damtraits.soils[order(damtoalphpop),,i][order(rep(Apopenvdat$TAnn,each=10)),])) 
		trtse <- as.vector(t(damtraitsse.soils[order(damtoalphpop),,i][order(rep(Apopenvdat$TAnn,each=10)),])) 
		plot(trtmns~rep(1:100,each=3),ylim=c(min(trtmns-trtse,na.rm=T),max(trtmns+trtse,na.rm=T)),
			ylab="",xlab="",xaxt="n",cex.axis=1.25,cex.lab=1,pch=NA )
		mtext(traitcols[i],side=1,line=0.5)
		lowy <- seq(from=0.5,to=100.5,by =10)
		for(j in 1:10){ 
			polygon(c(lowy[j],lowy[j]+10,lowy[j]+10,lowy[j]),rep(c(min(trtmns-trtse,na.rm=T),max(trtmns+trtse,na.rm=T)),each=2),border=NA,col=rbsoft2[order(Apopenvdat$TAnn)][j])
		}
		abline(v=c(lowy),lty=2,col=rgb(0,0,0,alpha=0.5))
		points(trtmns~rep(1:100,each=3), col=rep(c(1,rgb(0,1,0),rgb(1,0,1)), times=100),
			pch=rep(c(1,2,5),times=100), cex=.5)
		arrows(rep(1:100,each=3), trtmns-trtse, y1=trtmns+trtse,
			col=rep(c(rgb(0,0,0,alpha=0.5),rgb(0,1,0,alpha=0.5),rgb(1,0,1,alpha=0.5)),times=100),
			length=0 ,lwd=1 ) #in sel. skewers, mt is black, mc is green and mt is purple
		if(i==2){ legend(5,y=80,c("Biota15.0","Biota14.3","Biota13.0"),pch=c(1,2,5),col=c(rgb(0,0,0,alpha=0.5),rgb(0,1,0,alpha=0.5),rgb(1,0,1,alpha=0.5)),bty="n") }
	
		trtmnsp <- as.vector(t(poptraits.soils[order(Apopenvdat$TAnn),,i])) #mt first, them mc then tc
		trtsep <- as.vector(t(poptraitsse.soils[order(Apopenvdat$TAnn),,i])) 
		plot(trtmnsp~rep(1:10,each=3),ylim=c(min(trtmnsp)-max(trtsep),max(trtmnsp)+max(trtsep)),
			ylab="",xlab="",xaxt="n",cex.axis=1,cex.lab=1.25,pch=NA,xlim=c(0,11) )
		points(trtmnsp~rep(1:10,each=3),col=rep(c(rgb(0,0,0),rgb(0,1,0),rgb(1,0,1)), times=10),pch=rep(c(1,2,5),times=10),cex=.5 )
		arrows(rep(1:10,each=3),trtmnsp-trtsep,y1=trtmnsp+trtsep,col=rep(c(rgb(0,0,0,alpha=.5),rgb(0,1,0,alpha=.5),rgb(1,0,1,alpha=.5)), times=10),length=0 ,lwd=1 )
		lowy <- seq(from=0.5,to=10.5,by =1)
		for(j in 1:10){ 
			polygon(c(lowy[j],lowy[j]+1,lowy[j]+1,lowy[j]),rep(c(min(trtmnsp-trtsep,na.rm=T),max(trtmnsp+trtsep,na.rm=T)),each=2),border=NA,col=rbsoft2[order(Apopenvdat$TAnn)][j])
		}
		smns <- c(mt.mn[i],mc.mn[i],tc.mn[i])
		sses <- c(mt.se[i],mc.se[i],tc.se[i])
		plot(smns~c(1:3),ylim=c(min(smns)-max(sses),max(smns)+max(sses)),
		ylab="",xlab="",cex.axis=1,cex.lab=1.25,pch=NA ,xlim=c(0,4),xaxt="n")
		points(smns~c(1:3),col=rep(c(rgb(0,0,0),rgb(0,1,0),rgb(1,0,1)), times=10),pch=rep(c(1,2,5),times=10),cex=.5 )
		arrows(1:3,smns-sses,y1=smns+sses,col=rep(c(rgb(0,0,0,alpha=.5),rgb(0,1,0,alpha=.5),rgb(1,0,1,alpha=.5)), times=10),length=0 ,lwd=1 )

	}
dev.off() 
