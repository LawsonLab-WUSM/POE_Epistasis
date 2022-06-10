

rm(list = ls())

set.seed(0)

library(units)
library(sf)
library(spdep)
library(adegenet)
library(pegas)




getGenotype<-function(ImpScore,AddScore){
		
	if(ImpScore == 0 & AddScore == 1){
		return("L.L") ##LL
	}

	else if(ImpScore == 1 & AddScore == 0){
		return("L.S") ##LS
	}

	else if(ImpScore == -1 & AddScore == 0){
		return("S.L") ##SL
	}

	else if(ImpScore == 0 & AddScore == -1){
		return("S.S") ##SS
	}
	
	else {
		return("0.0")
	}
}



ImprintingScores<-read.delim("F16ImprintingScores.csv",sep=",",header=TRUE)
AdditiveScores<-read.delim("F16AdditiveScores.csv",sep=",",header=TRUE)

plinkIntermediate<-read.delim("plinkIntermediateF16Phenotypes.txt",sep="\t",header=FALSE)
colnames(plinkIntermediate)<-c("Family","Individual","PaternalID","MaternalID","Sex","Phenotype")


mapData<-read.delim("F16MapFile.map",header=FALSE)
colnames(mapData)<-c("Chr","marker","cM","Position.mm9")



tempImpScores<-ImprintingScores[,c(4,match(paste("X",plinkIntermediate$Individual,sep=""),colnames(ImprintingScores)))]
subetImpScores<-subset(tempImpScores,Marker.ID %in% mapData$marker)

tempAddScores<-AdditiveScores[,c(4,match(paste("X",plinkIntermediate$Individual,sep=""),colnames(AdditiveScores)))]
subetAddScores<-subset(tempAddScores,Marker.ID %in% mapData$marker)

subetImpScores[1]==subetAddScores[1]

genotypes<-subetImpScores

subetImpScores[,2:ncol(subetImpScores)]<-apply(subetImpScores[,2:ncol(subetImpScores)],2,function(N) round(N,digits=0))
subetAddScores[,2:ncol(subetAddScores)]<-apply(subetAddScores[,2:ncol(subetAddScores)],2,function(N) round(N,digits=0))

genotypes[2:ncol(subetImpScores)]<-data.frame(matrix( unlist(lapply(2:ncol(subetImpScores),function(I) lapply(1:nrow(subetImpScores), function(M) getGenotype(subetImpScores[M,I],subetAddScores[M,I]) ) )), byrow=FALSE, ncol=(ncol(subetImpScores)-1)))


GenotypeTrans<-t(genotypes[2:ncol(genotypes)])

colnames(GenotypeTrans)<-unlist(genotypes[1])


sum(paste("X",plinkIntermediate$Individual,sep="")!=rownames(GenotypeTrans))

pedData<-data.frame(plinkIntermediate, GenotypeTrans)

sum(mapData$marker!=colnames(GenotypeTrans))



splitGenotype<-function(Genotype){
	
	if(Genotype=="L.L"){return(c(1,1))}
	if(Genotype=="L.S"){return(c(1,2))}
	if(Genotype=="S.L"){return(c(2,1))}
	if(Genotype=="S.S"){return(c(2,2))}
	if(Genotype=="0.0"){return(c(NA,NA))}
}


TestLD2<-function(marker1, marker2){
	
	Marker.1.Column<-colnames(GenotypeTrans)==as.character(marker1)
	Genotypes1<-matrix(unlist(lapply(GenotypeTrans[,Marker.1.Column],splitGenotype)),byrow=TRUE,ncol=2)
	
	Marker.2.Column<-colnames(GenotypeTrans)==as.character(marker2)
	Genotypes2<-matrix(unlist(lapply(GenotypeTrans[,colnames(GenotypeTrans)==as.character(marker2)],splitGenotype)),byrow=TRUE,ncol=2)
	
	
	
	PairedDF<-data.frame(Genotypes1, Genotypes2)

	s <- apply(PairedDF, 1, anyNA)

	Loci.PairedDF <- alleles2loci(PairedDF[!s, ])
	LD2(Loci.PairedDF,locus=1:2)[2]
	
}


##### White adipose


Pairs.WAT<-read.delim("WAT_Marker_Pairs.txt",header=FALSE)
Pairs.WAT<-unique(Pairs.WAT)

#### All markers represented?
sum(!(unique(as.character(Pairs.WAT[,1])) %in% colnames(GenotypeTrans)))

sum(!(unique(as.character(Pairs.WAT[,2])) %in% colnames(GenotypeTrans)))


LD.Data.WAT<-data.frame(matrix(unlist(lapply(1:nrow(Pairs.WAT), function(R) TestLD2(Pairs.WAT[R,1],Pairs.WAT[R,2]))),byrow=TRUE,ncol=3))

colnames(LD.Data.WAT)<-c("T2","df","P-val")

LD.Data.WAT$fdr<-p.adjust(LD.Data.WAT$"P-val",method="fdr")

hist(LD.Data.WAT$"P-val")
hist(LD.Data.WAT$"fdr")

head(LD.Data.WAT)

LD.Data.WAT$Marker1<-Pairs.WAT[,1]
LD.Data.WAT$Marker2<-Pairs.WAT[,2]


head(LD.Data.WAT)

dim(LD.Data.WAT[LD.Data.WAT$fd<0.05,])[1]

dim(Pairs.WAT)[1]

dim(LD.Data.WAT[LD.Data.WAT$fd<0.05,])[1]/dim(Pairs.WAT)[1]



write.table(LD.Data.WAT,"Linkage_Statistics_WAT.txt",sep="\t", quote=FALSE,col.names=TRUE,row.names=FALSE)


#### Hypo

Pairs.Hypo<-read.delim("Hypo_Marker_Pairs.txt",header=FALSE)
Pairs.Hypo<-unique(Pairs.Hypo)

#### All markers represented?
sum(!(unique(as.character(Pairs.Hypo[,1])) %in% colnames(GenotypeTrans)))

sum(!(unique(as.character(Pairs.Hypo[,2])) %in% colnames(GenotypeTrans)))


LD.Data.Hypo<-data.frame(matrix(unlist(lapply(1:nrow(Pairs.Hypo), function(R) TestLD2(Pairs.Hypo[R,1],Pairs.Hypo[R,2]))),byrow=TRUE,ncol=3))

colnames(LD.Data.Hypo)<-c("T2","df","P-val")

LD.Data.Hypo$fdr<-p.adjust(LD.Data.Hypo$"P-val",method="fdr")

hist(LD.Data.Hypo$"P-val")
hist(LD.Data.Hypo$"fdr")

head(LD.Data.Hypo)

LD.Data.Hypo$Marker1<-Pairs.Hypo[,1]
LD.Data.Hypo$Marker2<-Pairs.Hypo[,2]


head(LD.Data.Hypo)

dim(LD.Data.Hypo[LD.Data.Hypo$fd<0.05,])[1]

dim(Pairs.Hypo)[1]

dim(LD.Data.Hypo[LD.Data.Hypo$fd<0.05,])[1]/dim(Pairs.Hypo)[1]



write.table(LD.Data.Hypo,"Linkage_Statistics_Hypo.txt",sep="\t", quote=FALSE,col.names=TRUE,row.names=FALSE)


#### Liver

Pairs.Liver<-read.delim("Liver_Marker_Pairs.txt",header=FALSE)
Pairs.Liver<-unique(Pairs.Liver)

#### All markers represented?
sum(!(unique(as.character(Pairs.Liver[,1])) %in% colnames(GenotypeTrans)))

sum(!(unique(as.character(Pairs.Liver[,2])) %in% colnames(GenotypeTrans)))



LD.Data.Liver<-data.frame(matrix(unlist(lapply(1:nrow(Pairs.Liver), function(R) TestLD2(Pairs.Liver[R,1],Pairs.Liver[R,2]))),byrow=TRUE,ncol=3))

colnames(LD.Data.Liver)<-c("T2","df","P-val")

LD.Data.Liver$fdr<-p.adjust(LD.Data.Liver$"P-val",method="fdr")

hist(LD.Data.Liver$"P-val")
hist(LD.Data.Liver$"fdr")

head(LD.Data.Liver)

LD.Data.Liver$Marker1<-Pairs.Liver[,1]
LD.Data.Liver$Marker2<-Pairs.Liver[,2]

TestLD2("wu-rs3671277","wu-rs3658906")

head(LD.Data.Liver)

dim(LD.Data.Liver[LD.Data.Liver$fd<0.05,])[1]

dim(Pairs.Liver)[1]

dim(LD.Data.Liver[LD.Data.Liver$fd<0.05,])[1]/dim(Pairs.Liver)[1]


write.table(LD.Data.Liver,"Linkage_Statistics_Liver.txt",sep="\t", quote=FALSE,col.names=TRUE,row.names=FALSE)


