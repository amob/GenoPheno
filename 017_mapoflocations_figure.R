library(geosphere)
library(raster)
library(ggmap)
library(rgdal)
library(dismo)

range01=function(x){
newnums <- sapply(1:length(x), function(z) ifelse(!is.na(x[z]), (x[z]-min(x,na.rm=T))/(max(x,na.rm=T)-min(x,na.rm=T)), NA ))
return(newnums)
}
opar <- par()


#get bioclim data
sites <- read.csv("~/Dropbox/AO-3 mine/population locations.csv",header=T,sep="\t")
#sort site file in meaningful way
sites$popabbr <- c("TZ","FP","MT","DA","MC","M","TC","TX","CL","CU")
sitesO <- sites[order(sites$popabbr),]
Ptpops <- data.frame(lon=sitesO$lon,lat=sitesO$lat)
SPtpops <- SpatialPoints(Ptpops,proj4string=CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"))
#read in ENV data!
files <- list.files(path ="~/bio_22_tif/", pattern=".tif", full.names = T)
bioclimvals <- matrix(nrow=nrow(Ptpops),ncol=length(files))
for(i in 1:length(files)){
	bioclimvals[,i]=extract(raster(files[i]),SPtpops) #bioclim crs is essentially crs="+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"
}
popbiodat <- data.frame(bioclimvals) #sorted by how are in sitesO, so alphabetical.
#columns are what the bioclim levels are(1 through 19): TAnn TmeanWarmQ TmeanCQ Pann PwetM PDM Pseas PwetQ PDQ PWarmQ PCQ TDiRange TIsother  TSD TmaxWarmM TminCM TannR TmeanWetQ TmeanDQ

MxCity <- SpatialPoints(data.frame(-99.1333,19.4328),proj4string=CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"))
#pts.ext <- extent(c(-106,-98,18,22))


mxalt <- crop(raster("~/maketeofigure/alt_22.tif"),extent(-116,-94.4,15,30)) #check if projection crs="+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"
mxsub <- crop(mxalt,extent(-100.25,-98.25,18.6,20))

bw <- colorRampPalette(c(rgb(0,0,0),rgb(1,1,1)))

par(opar)
pdf("~/Dropbox/AO-3 mine/PopLocs.pdf", width=6,height=4.85)
	plot(mxalt,col=bw(1000))
	mtext("Latitude",side=2,line=2)
	mtext("Longitude",side=1,line=2)
	points(MxCity, cex=2,pch='*',col=rgb(1,1,0))
	points(SPtpops, col = rgb(range01(popbiodat[,1]),0,1-(range01(popbiodat[,1]))), pch = 16,cex=0.5,zlim=c(0,5000))
	par(new = TRUE)
	par(fig = c(0.15, 0.55, 0.25, 0.65)) #c(x1, x2, y1, y2)  #my current x y sub extent spans 1/10 of the full xy extent, so these should be in equal units
	plot(mxsub,col=bw(1000),legend=FALSE,zlim=c(0,5000))
	par(fig = c(0.15, 0.55, 0.25, 0.65))
	points(MxCity, cex=2,pch='*',col=rgb(1,1,0))
	par(fig = c(0.15, 0.55, 0.25, 0.65))#par(fig = c(0.21, .59, 0.3, .5))
	points(SPtpops, col = rgb(range01(popbiodat[,1]),0,1-(range01(popbiodat[,1]))), pch = 16)
	points(SPtpops[6:8,], col = rgb(range01(popbiodat[,1])[6:8],0,1-range01(popbiodat[,1])[6:8]), pch = 1, cex=2)
dev.off()

#precip, temp, from worldclim by month for tile 22. 
files2 <- list.files(path ="~/Documents/prec_22_tif/", pattern=".tif", full.names = T)#order 1 10 11 12 2-9
precipmonth <- data.frame(matrix(nrow=10,ncol=12))
for(i in 1:12){
	x<-raster(files2[i])
	projection(x) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
	precipmonth[,i]<-extract(x,SPtpops)
}
colnames(precipmonth) <- c("Jan","Oct","Nov","Dec","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep")
precipmonth<-precipmonth[,c(1,5:12,2:4)]
#
files3 <- list.files(path ="~/Documents/tmean_22_tif/", pattern=".tif", full.names = T)#order 1 10 11 12 2-9
tmeanmonth <- data.frame(matrix(nrow=10,ncol=12))
for(i in 1:12){
	x<-raster(files3[i])
	projection(x) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
	tmeanmonth[,i]<-extract(x,SPtpops)
}
colnames(tmeanmonth) <- c("Jan","Oct","Nov","Dec","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep")
tmeanmonth<-tmeanmonth[,c(1,5:12,2:4)]

rb <-colorRampPalette(c(rgb(0,0,1),rgb(1,0,0)))
tempnames <- c("19.8","18.6","15.3","15.0","14.7","14.4","14.3","13.2","13.0","12.9")

pdf("~/Dropbox/AO-3 mine/PopMonthlyTempPrecip.pdf",width=5,height=5)#
	layout(matrix(1:4,nrow=2,byrow=T),widths=c(2.75,1))
	par(oma=c(2,0,1,1))
#	layout(matrix(1:5,ncol=5),widths=c(1.6,3,3,3,1.5))
	par(mar=c(1,4,0,1))
	matplot(t(precipmonth),type="l",col = rgb(range01(popbiodat[,1]),0,1-(range01(popbiodat[,1]))),lty=1,ylab="Total Monthly Precipitation mm",xaxt="n")
		abline(v=c(4,10),lty=2)#apprx germ and flowering months
	par(mar=c(1,6,3,0))
	image(1,seq(from=.5,to=9.5,by=1),t(matrix(rev(sort(popbiodat[,1])),nrow=10)),col=rb(50),zlim=c(129,198),ylim=c(0,10),ylab="",xlab="",bty="n",xaxt="n",yaxt="n")
		axis(side=2, at=seq(from=.5,to=9.5,by=1),labels=tempnames, las=2, cex.axis = 1)#
		mtext("Plant population MAT", side=2, line =3.75, cex =1)
	par(mar=c(1,4,0,1))
	matplot(t(tmeanmonth/10),type="l",col = rgb(range01(popbiodat[,1]),0,1-(range01(popbiodat[,1]))),lty=1,ylab="Montly Mean Temperature",xaxt="n")
		abline(v=c(4,10),lty=2)#apprx germ and flowering months
		axis(side=1,at=c(1:12),labels=c("J","F","M","A","M","J","J","A","S","O","N","D"))
dev.off()
#https://www.usda.gov/oce/weather/pubs/Other/MWCACP/Graphs/Mexico/MexSummerCornProd_0509.pdf
#assume planting = first possible germination, first harvest = first to seed = drydown
