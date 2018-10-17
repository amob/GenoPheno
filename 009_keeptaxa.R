names <- read.csv("~/GenoPheno/swarts_taxatokeep_firstpass.csv",stringsAsFactors=F)

missing <- read.csv("~/GenoPheno/missing_swarts2017.imiss",sep="\t",stringsAsFactors=F)

missingAOmex <- read.csv("~/GenoPheno/missing_AOMexGBS2.imiss",sep="\t",stringsAsFactors=F)

missinginset <-unlist(sapply(1:nrow(names), function(z) ifelse(names[z,1]%in%missing$INDV,missing$F_MISS[which(missing$INDV==names[z,1])],NA)))

namesout <- names[missinginset < .7 & !is.na(missinginset),1]
AOmexkeep <- missingAOmex$INDV[missingAOmex$F_MISS < .7]
# keep individuals from both sets that have missing data under 70%

write.table(namesout,"~/GenoPheno/SwartsKeepTaxa_indFilt70.txt",sep="\t",row.names=F,quote=F,col.names=F)

write.table(AOmexkeep, "~/GenoPheno/AOMexGBS2_indFilt70.txt",sep="\t",row.names=F,quote=F,col.names=F)
