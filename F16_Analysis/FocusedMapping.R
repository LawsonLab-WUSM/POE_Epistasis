
set.seed(0)
rm(list=ls())

library(lsr)
library(qvalue)
library(ggplot2)
library(cowplot)


#### Put File Names Below
Netmm9_BEDFile<-"WATNet_mm9.BED" ### Transcript Positions Bed file
NetAdditiveScoresFile<-"WATNetAdditiveScores.csv" ### F16 Additive Scores
NetImprintingScoresFile<-"WATNetImprintingScores.csv"	### F16 Imprinting Scores
PhenotypesFile<-"F16Phenotypes.csv"	### F16 Phenotypes
NetTableFile<-"WATNetTable.tsv" ### Relavent DE~ASE gene pairs AND relavent Phenotypes from QTL DE gene falls in
SNPBEDFile<-"F16_SNPS_within_3Mb_of_WATNet.BED" ### List of markers for each transcript

GenesInPOEQTL<-as.character(unique(read.table("genes_in_POE_QTL.BED")[,4]))

#### SNPBEDFile format: [Marker Chr] [Marker Start] [Marker Stop] [Transcript ID] [Gene Name]

#### NetTableFile format: [DE Gene] [ASE Gene] [pvalue] [F1 Trait Name] [F16 Trait Name]

##Import Data
{

Netmm9<-read.delim(Netmm9_BEDFile,header=FALSE,sep="\t")
colnames(Netmm9)<-c("Chr","Start","Stop","ENSEMBL","Gene")

### This is identifying the center point of each gene
NetMedian<-data.frame(matrix(unlist( lapply(1:length(unique(Netmm9$Gene)), function(n) c( mean(c(min(unlist(subset(Netmm9,Gene==unique(Netmm9$Gene)[n])[2:3])),max(unlist(subset(Netmm9,Gene==unique(Netmm9$Gene)[n])[2:3])))), as.character(unique(Netmm9$Gene)[n]) ) )), byrow=TRUE, ncol=2))

colnames(NetMedian)<-c("Locus","Gene")

NetAdditiveScores<-read.csv(NetAdditiveScoresFile)
NetAdditiveScores$LocusID<-apply(NetAdditiveScores[,c(2,5)],1,paste,collapse="_")
NetAdditiveScores$Marker.ID<-as.character(NetAdditiveScores$Marker.ID)

NetImprintingScores<-read.csv(NetImprintingScoresFile)
NetImprintingScores$LocusID<-apply(NetImprintingScores[,c(2,5)],1,paste,collapse="_")
NetImprintingScores$Marker.ID<-as.character(NetImprintingScores$Marker.ID)

F16Phenotypes<-read.csv(PhenotypesFile)

Net<-read.delim(NetTableFile,header=TRUE,sep="\t")

unique(Net[,4:5])


#Net <-subset(Net, Target.Gene %in% GenesInPOEQTL)

#nrow(Net)
#Net<-subset(Net,Trait %in% c("GlucWt1","GlucWt2","CHOL","TG","NecWeight","Insulin","AUC2","AUC1",'T0_1',"T0_2"))



GeneSNPs<-read.delim(SNPBEDFile,header=FALSE,sep="\t")
colnames(GeneSNPs)<-c("Chr","Start","Stop","ENSEMBL","Gene")

GeneSNPs$LocusID<-apply(GeneSNPs[1:2],1,paste,collapse="_")

###### SUPER IMPORTANT!!!!! Make sure R is treating phenotypes as numbers NOT factors
F16Phenotypes$Animal<-unlist(lapply(F16Phenotypes$Animal,function(Z) paste(c("X",Z),collapse="")))

}



### Transforming phenotype
{

plotTransform<-function(Phenotype,Original,Transformed){
	png(filename=paste(Phenotype,"_TransformationPlot.png",sep=""), width=500,height=500,type=c("cairo"))
	par(mfrow=c(2,2))
	
	hist(Original[,colnames(Original)==Phenotype],ylab="Frequency",xlab=Phenotype,main="Original")
	hist(Transformed[,colnames(Original)==Phenotype],ylab="Frequency",xlab=Phenotype,main="Transformed")
	
	qqnorm(Original[,colnames(Original)==Phenotype],pch=16,cex=0.5,main="Original")
	qqline(Original[,colnames(Original)==Phenotype],col="red")
	
	qqnorm(Transformed[,colnames(Original)==Phenotype],pch=16,cex=0.5,main="Tranformed")
	qqline(Transformed[,colnames(Original)==Phenotype],col="red")
	dev.off()
}



tF16Phenotypes<-F16Phenotypes

apply( tF16Phenotypes[,c(10:11,13:41)],2,function(p) shapiro.test(as.numeric(p[!is.na(p)]))[2]<0.05 )


colnames(tF16Phenotypes)


#tF16Phenotypes$Heart<-log10(F16Phenotypes$Heart) ### No transform!

#tF16Phenotypes$L.Kid<-log10(F16Phenotypes$L.Kid) ### No transform!
#tF16Phenotypes$L.Kid[!is.finite(tF16Phenotypes$L.Kid)]<-NA

#tF16Phenotypes$R.Kid<-log10(F16Phenotypes$R.Kid) ### No transform!
#tF16Phenotypes$R.Kid[!is.finite(tF16Phenotypes$R.Kid)]<-NA

tF16Phenotypes$Kidney<-log10(F16Phenotypes$Kidney) ### No transform!
tF16Phenotypes$Kidney[!is.finite(tF16Phenotypes$Kidney)]<-NA

tF16Phenotypes$Spleen<-log10(F16Phenotypes$Spleen)

tF16Phenotypes$Liver<-log10(F16Phenotypes$Liver)



tF16Phenotypes$GlucWt1<-log10(F16Phenotypes$GlucWt1)
#tF16Phenotypes$GlucWt2<-log10(F16Phenotypes$GlucWt2)

tF16Phenotypes$NecWeight<-log10(F16Phenotypes$NecWeight)
#tF16Phenotypes$NecWeight<-scale(F16Phenotypes$NecWeight,scale=TRUE,center=TRUE)[,1]

tF16Phenotypes$Go.dal<-log10(F16Phenotypes$Go.dal)
tF16Phenotypes$Kid_Fat<-log10(F16Phenotypes$Kid_Fat)

#tF16Phenotypes$Mesen<-log10(F16Phenotypes$Mesen) ### No transform!

#tF16Phenotypes$Ing<-log10(F16Phenotypes$Ing) ### No transform!
#tF16Phenotypes$Ing[!is.finite(tF16Phenotypes$Ing)]<-NA

#tF16Phenotypes$Tot_Fat<-log2(F16Phenotypes$Tot_Fat) ### No transform!
#tF16Phenotypes$Tot_Fat<-scale(F16Phenotypes$Tot_Fat,scale=TRUE,center=TRUE) ### No transform!



#tF16Phenotypes$T0_1<-log10(F16Phenotypes$T0_1) ### No transform!
tF16Phenotypes$T15_1<-log10(F16Phenotypes$T15_1)
tF16Phenotypes$T30_1<-log10(F16Phenotypes$T30_1)
tF16Phenotypes$T60_1<-log10(F16Phenotypes$T60_1)
tF16Phenotypes$T120_1<-log10(F16Phenotypes$T120_1)

tF16Phenotypes$AUC1<-log10(F16Phenotypes$AUC1)
tF16Phenotypes$AUC1[!is.finite(tF16Phenotypes$AUC1)]<-NA

tF16Phenotypes$AUC2<-log10(F16Phenotypes$AUC2)
tF16Phenotypes$AUC2[!is.finite(tF16Phenotypes$AUC2)]<-NA

#tF16Phenotypes$T0_2<-log10(F16Phenotypes$T0_2) ### No transform!
tF16Phenotypes$T15_2<-log10(F16Phenotypes$T15_2)
tF16Phenotypes$T30_2<-log10(F16Phenotypes$T30_2)
tF16Phenotypes$T60_2<-log10(F16Phenotypes$T60_2)
tF16Phenotypes$T120_2<-log10(F16Phenotypes$T120_2)

tF16Phenotypes$GLC<-log10(F16Phenotypes$GLC)
tF16Phenotypes$Insulin<-log10(F16Phenotypes$Insulin)

#tF16Phenotypes$FFA<-log2(F16Phenotypes$FFA) ### No transform!
#tF16Phenotypes$FFA[!is.finite(tF16Phenotypes$FFA)]<-NA

#tF16Phenotypes$CHOL<-log10(F16Phenotypes$CHOL) ### No transform!

tF16Phenotypes$TG<-log10(F16Phenotypes$TG)
tF16Phenotypes$TG[!is.finite(tF16Phenotypes$TG)]<-NA



apply( F16Phenotypes[,c(10:11,13:41)],2,function(p) shapiro.test(as.numeric(p[!is.na(p)]))[2]<0.05 )

unlist(apply( F16Phenotypes[,c(10:11,13:41)],2,function(p) shapiro.test(as.numeric(p[!is.na(p)]))[2] )) - unlist(apply( tF16Phenotypes[,c(10:11,13:41)],2,function(p) shapiro.test(as.numeric(p[!is.na(p)]))[2] ))

apply( tF16Phenotypes[,c(10:11,13:41)],2,function(p) shapiro.test(as.numeric(p[!is.na(p)]))[2]<0.05 )


plotTransform("FFA",F16Phenotypes,tF16Phenotypes)
plotTransform("CHOL",F16Phenotypes,tF16Phenotypes)
plotTransform("TG",F16Phenotypes,tF16Phenotypes)


plotTransform("GlucWt1",F16Phenotypes,tF16Phenotypes)
plotTransform("GlucWt2",F16Phenotypes,tF16Phenotypes)

plotTransform("NecWeight",F16Phenotypes,tF16Phenotypes)
plotTransform("Go.dal",F16Phenotypes,tF16Phenotypes)
plotTransform("Kid_Fat",F16Phenotypes,tF16Phenotypes)
plotTransform("Mesen",F16Phenotypes,tF16Phenotypes)
plotTransform("Ing",F16Phenotypes,tF16Phenotypes)
plotTransform("Tot_Fat",F16Phenotypes,tF16Phenotypes)


plotTransform("Heart",F16Phenotypes,tF16Phenotypes)
plotTransform("L.Kid",F16Phenotypes,tF16Phenotypes)
plotTransform("R.Kid",F16Phenotypes,tF16Phenotypes)
plotTransform("Kidney",F16Phenotypes,tF16Phenotypes)
plotTransform("Spleen",F16Phenotypes,tF16Phenotypes)
plotTransform("Liver",F16Phenotypes,tF16Phenotypes)


plotTransform("Insulin",F16Phenotypes,tF16Phenotypes)
plotTransform("GLC",F16Phenotypes,tF16Phenotypes)

plotTransform("T0_1",F16Phenotypes,tF16Phenotypes)
plotTransform("T15_1",F16Phenotypes,tF16Phenotypes)
plotTransform("T30_1",F16Phenotypes,tF16Phenotypes)
plotTransform("T60_1",F16Phenotypes,tF16Phenotypes)
plotTransform("T120_1",F16Phenotypes,tF16Phenotypes)

plotTransform("T0_2",F16Phenotypes,tF16Phenotypes)
plotTransform("T15_2",F16Phenotypes,tF16Phenotypes)
plotTransform("T30_2",F16Phenotypes,tF16Phenotypes)
plotTransform("T60_2",F16Phenotypes,tF16Phenotypes)
plotTransform("T120_2",F16Phenotypes,tF16Phenotypes)

plotTransform("AUC1",F16Phenotypes,tF16Phenotypes)
plotTransform("AUC2",F16Phenotypes,tF16Phenotypes)

}


### Covariate Screen
{

head(tF16Phenotypes)
sF16Phenotypes<-tF16Phenotypes

sF16Phenotypes$GlucWt1[as.numeric(rownames(data.frame(residuals(aov(GlucWt1~Sex*Diet*+Litter_Size+Family,sF16Phenotypes)))))]<-residuals(aov(GlucWt1~Sex*Diet*+Litter_Size+Family,sF16Phenotypes))

sF16Phenotypes$GlucWt2[as.numeric(rownames(data.frame(residuals(aov(GlucWt2~Sex*Diet*+Litter_Size+Family,sF16Phenotypes)))))]<-residuals(aov(GlucWt2~Sex*Diet*+Litter_Size+Family,sF16Phenotypes))

sF16Phenotypes$T0_1[as.numeric(rownames(data.frame(residuals(aov(T0_1~Sex*Diet*+Litter_Size+Family,sF16Phenotypes)))))]<-residuals(aov(T0_1~Sex*Diet*+Litter_Size+Family,sF16Phenotypes))

sF16Phenotypes$T15_1[as.numeric(rownames(data.frame(residuals(aov(T15_1~Sex*Diet*+Litter_Size+Family,sF16Phenotypes)))))]<-residuals(aov(T15_1~Sex*Diet*+Litter_Size+Family,sF16Phenotypes))

sF16Phenotypes$T30_1[as.numeric(rownames(data.frame(residuals(aov(T30_1~Sex*Diet*+Litter_Size+Family,sF16Phenotypes)))))]<-residuals(aov(T30_1~Sex*Diet*+Litter_Size+Family,sF16Phenotypes))

sF16Phenotypes$T60_1[as.numeric(rownames(data.frame(residuals(aov(T60_1~Sex*Diet*+Litter_Size+Family,sF16Phenotypes)))))]<-residuals(aov(T60_1~Sex*Diet*+Litter_Size+Family,sF16Phenotypes))

sF16Phenotypes$T120_1[as.numeric(rownames(data.frame(residuals(aov(T120_1~Sex*Diet*+Litter_Size+Family,sF16Phenotypes)))))]<-residuals(aov(T120_1~Sex*Diet*+Litter_Size+Family,sF16Phenotypes))

sF16Phenotypes$AUC1[as.numeric(rownames(data.frame(residuals(aov(AUC1~Sex*Diet*+Litter_Size+Family,sF16Phenotypes)))))]<-residuals(aov(AUC1~Sex*Diet*+Litter_Size+Family,sF16Phenotypes))


sF16Phenotypes$T0_2[as.numeric(rownames(data.frame(residuals(aov(T0_2~Sex*Diet*+Litter_Size+Family,sF16Phenotypes)))))]<-residuals(aov(T0_2~Sex*Diet*+Litter_Size+Family,sF16Phenotypes))

sF16Phenotypes$T15_2[as.numeric(rownames(data.frame(residuals(aov(T15_2~Sex*Diet*+Litter_Size+Family,sF16Phenotypes)))))]<-residuals(aov(T15_2~Sex*Diet*+Litter_Size+Family,sF16Phenotypes))

sF16Phenotypes$T30_2[as.numeric(rownames(data.frame(residuals(aov(T30_2~Sex*Diet*+Litter_Size+Family,sF16Phenotypes)))))]<-residuals(aov(T30_2~Sex*Diet*+Litter_Size+Family,sF16Phenotypes))

sF16Phenotypes$T60_2[as.numeric(rownames(data.frame(residuals(aov(T60_2~Sex*Diet*+Litter_Size+Family,sF16Phenotypes)))))]<-residuals(aov(T60_2~Sex*Diet*+Litter_Size+Family,sF16Phenotypes))

sF16Phenotypes$T120_2[as.numeric(rownames(data.frame(residuals(aov(T120_2~Sex*Diet*+Litter_Size+Family,sF16Phenotypes)))))]<-residuals(aov(T120_2~Sex*Diet*+Litter_Size+Family,sF16Phenotypes))

sF16Phenotypes$AUC2[as.numeric(rownames(data.frame(residuals(aov(AUC2~Sex*Diet*+Litter_Size+Family,sF16Phenotypes)))))]<-residuals(aov(AUC2~Sex*Diet*+Litter_Size+Family,sF16Phenotypes))


sF16Phenotypes$NecWeight[as.numeric(rownames(data.frame(residuals(aov(NecWeight~Sex*Diet*+Litter_Size+Family,sF16Phenotypes)))))]<-residuals(aov(NecWeight~Sex*Diet*+Litter_Size+Family,sF16Phenotypes))

sF16Phenotypes$Heart[as.numeric(rownames(data.frame(residuals(aov(Heart~Sex*Diet*+Litter_Size+Family,sF16Phenotypes)))))]<-residuals(aov(Heart~Sex*Diet*+Litter_Size+Family,sF16Phenotypes))

sF16Phenotypes$L.Kid[as.numeric(rownames(data.frame(residuals(aov(L.Kid~Sex*Diet*+Litter_Size+Family,sF16Phenotypes)))))]<-residuals(aov(L.Kid~Sex*Diet*+Litter_Size+Family,sF16Phenotypes))

sF16Phenotypes$R.Kid[as.numeric(rownames(data.frame(residuals(aov(R.Kid~Sex*Diet*+Litter_Size+Family,sF16Phenotypes)))))]<-residuals(aov(R.Kid~Sex*Diet*+Litter_Size+Family,sF16Phenotypes))

sF16Phenotypes$Kidney[as.numeric(rownames(data.frame(residuals(aov(Kidney~Sex*Diet*+Litter_Size+Family,sF16Phenotypes)))))]<-residuals(aov(Kidney~Sex*Diet*+Litter_Size+Family,sF16Phenotypes))

sF16Phenotypes$Spleen[as.numeric(rownames(data.frame(residuals(aov(Spleen~Sex*Diet*+Litter_Size+Family,sF16Phenotypes)))))]<-residuals(aov(Spleen~Sex*Diet*+Litter_Size+Family,sF16Phenotypes))

sF16Phenotypes$Liver[as.numeric(rownames(data.frame(residuals(aov(Liver~Sex*Diet*+Litter_Size+Family,sF16Phenotypes)))))]<-residuals(aov(Liver~Sex*Diet*+Litter_Size+Family,sF16Phenotypes))

sF16Phenotypes$Go.dal[as.numeric(rownames(data.frame(residuals(aov(Go.dal~Sex*Diet*+Litter_Size+Family,sF16Phenotypes)))))]<-residuals(aov(Go.dal~Sex*Diet*+Litter_Size+Family,sF16Phenotypes))

sF16Phenotypes$Kid_Fat[as.numeric(rownames(data.frame(residuals(aov(Kid_Fat~Sex*Diet*+Litter_Size+Family,sF16Phenotypes)))))]<-residuals(aov(Kid_Fat~Sex*Diet*+Litter_Size+Family,sF16Phenotypes))

sF16Phenotypes$Mesen[as.numeric(rownames(data.frame(residuals(aov(Mesen~Sex*Diet*+Litter_Size+Family,sF16Phenotypes)))))]<-residuals(aov(Mesen~Sex*Diet*+Litter_Size+Family,sF16Phenotypes))

sF16Phenotypes$Ing[as.numeric(rownames(data.frame(residuals(aov(Ing~Sex*Diet*+Litter_Size+Family,sF16Phenotypes)))))]<-residuals(aov(Ing~Sex*Diet*+Litter_Size+Family,sF16Phenotypes))

sF16Phenotypes$FFA[as.numeric(rownames(data.frame(residuals(aov(FFA~Sex*Diet*+Litter_Size+Family,sF16Phenotypes)))))]<-residuals(aov(FFA~Sex*Diet*+Litter_Size+Family,sF16Phenotypes))

sF16Phenotypes$GLC[as.numeric(rownames(data.frame(residuals(aov(GLC~Sex*Diet*+Litter_Size+Family,sF16Phenotypes)))))]<-residuals(aov(GLC~Sex*Diet*+Litter_Size+Family,sF16Phenotypes))

sF16Phenotypes$CHOL[as.numeric(rownames(data.frame(residuals(aov(CHOL~Sex*Diet*+Litter_Size+Family,sF16Phenotypes)))))]<-residuals(aov(CHOL~Sex*Diet*+Litter_Size+Family,sF16Phenotypes))

sF16Phenotypes$TG[as.numeric(rownames(data.frame(residuals(aov(TG~Sex*Diet*+Litter_Size+Family,sF16Phenotypes)))))]<-residuals(aov(TG~Sex*Diet*+Litter_Size+Family,sF16Phenotypes))

sF16Phenotypes$Insulin[as.numeric(rownames(data.frame(residuals(aov(Insulin~Sex*Diet*+Litter_Size+Family,sF16Phenotypes)))))]<-residuals(aov(Insulin~Sex*Diet*+Litter_Size+Family,sF16Phenotypes))

sF16Phenotypes$Tot_Fat[as.numeric(rownames(data.frame(residuals(aov(Tot_Fat~Sex*Diet*+Litter_Size+Family,sF16Phenotypes)))))]<-residuals(aov(Tot_Fat~Sex*Diet*+Litter_Size+Family,sF16Phenotypes))

unique(apply(sF16Phenotypes[7:8],1,function(S) paste(S,collapse=".")))

cohortColors<-c("darkRed","red","black","darkBlue","blue")[as.numeric(as.factor(apply(sF16Phenotypes[7:8],1,function(S) paste(S,collapse="."))))]


png(filename=paste("GlucWt1","_ResidualPlot.png",sep=""), width=500,height=500,type=c("cairo"))
plot(sF16Phenotypes$GlucWt1, F16Phenotypes$GlucWt1,pch=16,cex=0.5,col=cohortColors,ylab="Raw Data",xlab="Residualized Tranformed Data",main="GlucWt1")
dev.off()

png(filename=paste("GlucWt2","_ResidualPlot.png",sep=""), width=500,height=500,type=c("cairo"))
plot(sF16Phenotypes$GlucWt2, F16Phenotypes$GlucWt2,pch=16,cex=0.5,col=cohortColors,ylab="Raw Data",xlab="Residualized Tranformed Data",main="GlucWt2")
dev.off()

png(filename=paste("T0_1","_ResidualPlot.png",sep=""), width=500,height=500,type=c("cairo"))
plot(sF16Phenotypes$T0_1, F16Phenotypes$T0_1,pch=16,cex=0.5,col=cohortColors,ylab="Raw Data",xlab="Residualized Tranformed Data",main="T0_1")
dev.off()

png(filename=paste("T15_1","_ResidualPlot.png",sep=""), width=500,height=500,type=c("cairo"))
plot(sF16Phenotypes$T15_1, F16Phenotypes$T15_1,pch=16,cex=0.5,col=cohortColors,ylab="Raw Data",xlab="Residualized Tranformed Data",main="T15_1")
dev.off()

png(filename=paste("T30_1","_ResidualPlot.png",sep=""), width=500,height=500,type=c("cairo"))
plot(sF16Phenotypes$T30_1, F16Phenotypes$T30_1,pch=16,cex=0.5,col=cohortColors,ylab="Raw Data",xlab="Residualized Tranformed Data",main="T30_1")
dev.off()

png(filename=paste("T60_1","_ResidualPlot.png",sep=""), width=500,height=500,type=c("cairo"))
plot(sF16Phenotypes$T60_1, F16Phenotypes$T60_1,pch=16,cex=0.5,col=cohortColors,ylab="Raw Data",xlab="Residualized Tranformed Data",main="T60_1")
dev.off()

png(filename=paste("T120_1","_ResidualPlot.png",sep=""), width=500,height=500,type=c("cairo"))
plot(sF16Phenotypes$T120_1, F16Phenotypes$T120_1,pch=16,cex=0.5,col=cohortColors,ylab="Raw Data",xlab="Residualized Tranformed Data",main="T120_1")
dev.off()

png(filename=paste("AUC1","_ResidualPlot.png",sep=""), width=500,height=500,type=c("cairo"))
plot(sF16Phenotypes$AUC1, F16Phenotypes$AUC1,pch=16,cex=0.5,col=cohortColors,ylab="Raw Data",xlab="Residualized Tranformed Data",main="AUC1")
dev.off()

png(filename=paste("T0_2","_ResidualPlot.png",sep=""), width=500,height=500,type=c("cairo"))
plot(sF16Phenotypes$T0_2, F16Phenotypes$T0_2,pch=16,cex=0.5,col=cohortColors,ylab="Raw Data",xlab="Residualized Tranformed Data",main="T0_2")
dev.off()

png(filename=paste("T15_2","_ResidualPlot.png",sep=""), width=500,height=500,type=c("cairo"))
plot(sF16Phenotypes$T15_2, F16Phenotypes$T15_2,pch=16,cex=0.5,col=cohortColors,ylab="Raw Data",xlab="Residualized Tranformed Data",main="T15_2")
dev.off()

png(filename=paste("T30_2","_ResidualPlot.png",sep=""), width=500,height=500,type=c("cairo"))
plot(sF16Phenotypes$T30_2, F16Phenotypes$T30_2,pch=16,cex=0.5,col=cohortColors,ylab="Raw Data",xlab="Residualized Tranformed Data",main="T30_2")
dev.off()

png(filename=paste("T60_2","_ResidualPlot.png",sep=""), width=500,height=500,type=c("cairo"))
plot(sF16Phenotypes$T60_2, F16Phenotypes$T60_2,pch=16,cex=0.5,col=cohortColors,ylab="Raw Data",xlab="Residualized Tranformed Data",main="T60_2")
dev.off()

png(filename=paste("T120_2","_ResidualPlot.png",sep=""), width=500,height=500,type=c("cairo"))
plot(sF16Phenotypes$T120_2, F16Phenotypes$T120_2,pch=16,cex=0.5,col=cohortColors,ylab="Raw Data",xlab="Residualized Tranformed Data",main="T120_2")
dev.off()

png(filename=paste("AUC2","_ResidualPlot.png",sep=""), width=500,height=500,type=c("cairo"))
plot(sF16Phenotypes$AUC2, F16Phenotypes$AUC2,pch=16,cex=0.5,col=cohortColors,ylab="Raw Data",xlab="Residualized Tranformed Data",main="AUC2")
dev.off()

png(filename=paste("NecWeight","_ResidualPlot.png",sep=""), width=500,height=500,type=c("cairo"))
plot(sF16Phenotypes$NecWeight, F16Phenotypes$NecWeight,pch=16,cex=0.5,col=cohortColors,ylab="Raw Data",xlab="Residualized Tranformed Data",main="NecWeight")
dev.off()

png(filename=paste("Heart","_ResidualPlot.png",sep=""), width=500,height=500,type=c("cairo"))
plot(sF16Phenotypes$Heart, F16Phenotypes$Heart,pch=16,cex=0.5,col=cohortColors,ylab="Raw Data",xlab="Residualized Tranformed Data",main="Heart")
dev.off()

png(filename=paste("L.Kid","_ResidualPlot.png",sep=""), width=500,height=500,type=c("cairo"))
plot(sF16Phenotypes$L.Kid, F16Phenotypes$L.Kid,pch=16,cex=0.5,col=cohortColors,ylab="Raw Data",xlab="Residualized Tranformed Data",main="L.Kid")
dev.off()

png(filename=paste("R.Kid","_ResidualPlot.png",sep=""), width=500,height=500,type=c("cairo"))
plot(sF16Phenotypes$R.Kid, F16Phenotypes$R.Kid,pch=16,cex=0.5,col=cohortColors,ylab="Raw Data",xlab="Residualized Tranformed Data",main="R.Kid")
dev.off()

png(filename=paste("Kidney","_ResidualPlot.png",sep=""), width=500,height=500,type=c("cairo"))
plot(sF16Phenotypes$Kidney, F16Phenotypes$Kidney,pch=16,cex=0.5,col=cohortColors,ylab="Raw Data",xlab="Residualized Tranformed Data",main="Kidney")
dev.off()

png(filename=paste("Spleen","_ResidualPlot.png",sep=""), width=500,height=500,type=c("cairo"))
plot(sF16Phenotypes$Spleen, F16Phenotypes$Spleen,pch=16,cex=0.5,col=cohortColors,ylab="Raw Data",xlab="Residualized Tranformed Data",main="Spleem")
dev.off()

png(filename=paste("Liver","_ResidualPlot.png",sep=""), width=500,height=500,type=c("cairo"))
plot(sF16Phenotypes$Liver, F16Phenotypes$Liver,pch=16,cex=0.5,col=cohortColors,ylab="Raw Data",xlab="Residualized Tranformed Data",main="Liver")
dev.off()

png(filename=paste("Go.dal","_ResidualPlot.png",sep=""), width=500,height=500,type=c("cairo"))
plot(sF16Phenotypes$Go.dal, F16Phenotypes$Go.dal,pch=16,cex=0.5,col=cohortColors,ylab="Raw Data",xlab="Residualized Tranformed Data",main="Go.dal")
dev.off()

png(filename=paste("Kid_Fat","_ResidualPlot.png",sep=""), width=500,height=500,type=c("cairo"))
plot(sF16Phenotypes$Kid_Fat, F16Phenotypes$Kid_Fat,pch=16,cex=0.5,col=cohortColors,ylab="Raw Data",xlab="Residualized Tranformed Data",main="Kid_Fat")
dev.off()

png(filename=paste("Mesen","_ResidualPlot.png",sep=""), width=500,height=500,type=c("cairo"))
plot(sF16Phenotypes$Mesen, F16Phenotypes$Mesen,pch=16,cex=0.5,col=cohortColors,ylab="Raw Data",xlab="Residualized Tranformed Data",main="Mesen")
dev.off()

png(filename=paste("Ing","_ResidualPlot.png",sep=""), width=500,height=500,type=c("cairo"))
plot(sF16Phenotypes$Ing, F16Phenotypes$Ing,pch=16,cex=0.5,col=cohortColors,ylab="Raw Data",xlab="Residualized Tranformed Data",main="Ing")
dev.off()

png(filename=paste("FFA","_ResidualPlot.png",sep=""), width=500,height=500,type=c("cairo"))
plot(sF16Phenotypes$FFA, F16Phenotypes$FFA,pch=16,cex=0.5,col=cohortColors,ylab="Raw Data",xlab="Residualized Tranformed Data",main="FFA")
dev.off()

png(filename=paste("GLC","_ResidualPlot.png",sep=""), width=500,height=500,type=c("cairo"))
plot(sF16Phenotypes$GLC, F16Phenotypes$GLC,pch=16,cex=0.5,col=cohortColors,ylab="Raw Data",xlab="Residualized Tranformed Data",main="GLC")
dev.off()

png(filename=paste("CHOL","_ResidualPlot.png",sep=""), width=500,height=500,type=c("cairo"))
plot(sF16Phenotypes$CHOL, F16Phenotypes$CHOL,pch=16,cex=0.5,col=cohortColors,ylab="Raw Data",xlab="Residualized Tranformed Data",main="CHOL")
dev.off()

png(filename=paste("TG","_ResidualPlot.png",sep=""), width=500,height=500,type=c("cairo"))
plot(sF16Phenotypes$TG, F16Phenotypes$TG,pch=16,cex=0.5,col=cohortColors,ylab="Raw Data",xlab="Residualized Tranformed Data",main="TG")
dev.off()

png(filename=paste("Insulin","_ResidualPlot.png",sep=""), width=500,height=500,type=c("cairo"))
plot(sF16Phenotypes$Insulin, F16Phenotypes$Insulin,pch=16,cex=0.5,col=cohortColors,ylab="Raw Data",xlab="Residualized Tranformed Data",main="Insulin")
dev.off()

png(filename=paste("Tot_Fat","_ResidualPlot.png",sep=""), width=500,height=500,type=c("cairo"))
plot(sF16Phenotypes$Tot_Fat, F16Phenotypes$Tot_Fat,pch=16,cex=0.5,col=cohortColors,ylab="Raw Data",xlab="Residualized Tranformed Data",main="Tot_Fat")
dev.off()


}










plotPermFDRQC<-function(PVALUES,QVALUES,QCUTOFFVALUE,NULL.PVALUES){

	
DF<-data.frame(PVALUES,QVALUES)




QvP<-ggplot(DF,aes(x= PVALUES, y= QVALUES))+
geom_line(size=1)+
xlab("p-values")+
ylab("q-values")+
geom_hline(yintercept=QCUTOFFVALUE,color="red")

DF2<-data.frame(
QvalueCutoff=seq(0,1,0.01),
SigTests=unlist(lapply(seq(0,1,0.01),function(CUTOFF) sum(DF$QVALUES<=CUTOFF,na.rm=TRUE) )),
FalsePositives=unlist(lapply(seq(0,1,0.01),function(CUTOFF) sum(DF$QVALUES<=CUTOFF,na.rm=TRUE) ))*seq(0,1,0.01)
)


# BREAKS<-seq(min(DF2$SigTests),max(DF2$SigTests),round(abs(max(DF2$SigTests)-min(DF2$SigTests))/10,0))

BREAKS<-NULL
INIT<-min(DF2$SigTests,na.rm=TRUE)+0.5
FIN<-max(DF2$SigTests,na.rm=TRUE)


while(INIT<=FIN){
	BREAKS <-c(BREAKS,round(INIT,0))
	INIT<-INIT*2
}


SigvQ<-ggplot(DF2,aes(x=QvalueCutoff, y=SigTests+0.1))+
geom_line(size=1)+
scale_y_continuous(name="Significant tests",trans=scales::pseudo_log_trans(base=2), breaks=c(BREAKS) )+
scale_x_continuous(name="q-value cutoff",breaks=seq(0,1,0.1) )+
geom_vline(xintercept=QCUTOFFVALUE,color="red")






CORRESPONDINGROW<-c(1:length(seq(0,1,0.01)))[seq(0,1,0.01)==QCUTOFFVALUE]

SigvFalse<-ggplot(DF2,aes(x=SigTests,y=FalsePositives))+
geom_line(size=1)+
scale_y_continuous(name="Expected false positives",trans=scales::pseudo_log_trans(base=2), breaks=c(BREAKS) )+
scale_x_continuous(name="Significant tests",trans=scales::pseudo_log_trans(base=2), breaks=c(BREAKS))+
geom_vline(xintercept=DF2[CORRESPONDINGROW,2],color="red")+
theme(axis.text.x = element_text(angle = 45, hjust = 1))





COMB.DF<-data.frame(
pvalue=c(NULL.PVALUES, PVALUES),
Group=c(rep("Null",length(NULL.PVALUES)),rep("Real",length(PVALUES)))
)


# BINWIDTH<-max(DF$PVALUES[DF$QVALUES<=QCUTOFFVALUE],na.rm=TRUE)*10
BINWIDTH<-QCUTOFFVALUE/10

print(BINWIDTH)

COMPARE.DISTS<-ggplot(COMB.DF,aes(pvalue))+
geom_histogram(data=subset(COMB.DF,Group=="Real" & !is.na(pvalue) ),aes(y=..count../sum(..count..)),fill="blue",alpha=0.6,binwidth= BINWIDTH)+
geom_histogram(data=subset(COMB.DF,Group=="Null" & !is.na(pvalue)),aes(y=..count../sum(..count..)),fill="orange",alpha=0.6,binwidth= BINWIDTH)+
scale_x_continuous(expand = c(0, 0))+
scale_y_continuous(expand = c(0, 0))+
ylab("P(pvalue)")


COMBINED<-plot_grid(COMPARE.DISTS, QvP, SigvQ, SigvFalse, ncol=2)

return(COMBINED)

	
}










Bar_Dir<-"Epistasis_Barplots"
dir.create(Bar_Dir)



### Epistasis Test Functions
{

AddAndImpTest<-function(locus1,locus2,ImprintingScores,AdditiveScores,Phenotypes,phenOfInterest,DEgene,ASEgene,RANDOM=FALSE){
	
	### Pulling matched phenotypes for individuals from Imprinting scores object
	matchedPhenotypes<-Phenotypes[match(colnames(ImprintingScores[locus1,8:ncol(ImprintingScores)]),Phenotypes$Animal),]
	
	### Pulling matched columns for phenotypes of interest
	matchedPhenIndex<-match(phenOfInterest,colnames(matchedPhenotypes))
	
	### Construct data frame from marker genotypes and matched phenotypes
	datMat<-cbind(t(data.frame(ImprintingScores[c(locus1,locus2),8:ncol(ImprintingScores)])), matchedPhenotypes[,matchedPhenIndex])
	
	
	colnames(datMat)<-c("L1_Imp", "L2_Imp","phen")
	datMat<-data.frame(apply(datMat,2,as.numeric))
	
	if(RANDOM == TRUE){
		datMat$phen <-sample(datMat$phen,nrow(datMat))
	}
	
	
	
	### Through out any animals with no phenotype
	datMat<-datMat[!is.na(datMat$phen),]
	
	### Run ANOVA model and return F-value's 
	c( summary(aov(phen~.*.,data.frame(datMat)))[[1]][["Pr(>F)"]][1:3],DEgene,ASEgene,phenOfInterest,ImprintingScores[locus1,4],ImprintingScores[locus2,4] )
	
}





### This runs the test for All combinations of a gene Pairs markers for the relavent phenotype
genePairEpiTest<-function(gDE,gImp,SNPS,ImprintingScores,AdditiveScores,Phenotypes,phenOfInterest,RANDOM=FALSE){

	### Pulls all loci for DE Gene
	gDELoci<-match(unique(subset(SNPS,Gene==gDE)$LocusID), ImprintingScores$LocusID)
	
	### Pulls all loci for ASE Gene
	gImpLoci<-match(unique(subset(SNPS,Gene== gImp)$LocusID), ImprintingScores$LocusID)
	
	### Generates matrix of All Loci Pairs
	LociMat<-matrix(unlist(lapply(gDELoci, function(A) lapply(gImpLoci,function(B) return(c(A,B))) )),byrow=TRUE,ncol=2)
	
	### Runs Test on all pairs in matrix. Returns list of F-statistics
	unlist((apply(LociMat,1,function(X) AddAndImpTest(X[1], X[2], ImprintingScores, AdditiveScores, Phenotypes, phenOfInterest, gDE, gImp, RANDOM))))

}







PlotImpTest<-function(locus1,locus2,ImprintingScores,AdditiveScores,Phenotypes,phenOfInterest,DEgene,ASEgene,BARDIR){
	
	### Pulling matched phenotypes for individuals from Imprinting scores object
	matchedPhenotypes<-Phenotypes[match(colnames(ImprintingScores[locus1,8:ncol(ImprintingScores)]),Phenotypes$Animal),]
	
	### Pulling matched columns for phenotypes of interest
	matchedPhenIndex<-match(phenOfInterest,colnames(matchedPhenotypes))
	
	### Construct data frame from marker genotypes and matched phenotypes
	datMat<-cbind(t(data.frame(ImprintingScores[c(locus1,locus2),8:ncol(ImprintingScores)])), matchedPhenotypes[,matchedPhenIndex])
	
	
	
	
	
	
	
	colnames(datMat)<-c("L1_Imp", "L2_Imp","phen")
	datMat<-data.frame(apply(datMat,2,as.numeric))
	
	### Through out any animals with no phenotype
	datMat<-datMat[!is.na(datMat$phen),]
	

	DFMEANS<-data.frame(
	L1_Imp=c(1,-1,1,-1),	
	L2_Imp=c(1,1,-1,-1),
	Means=c(
	mean(subset(datMat, L1_Imp == 1 & L2_Imp == 1 )$phen),
	mean(subset(datMat, L1_Imp == -1 & L2_Imp == 1 )$phen),
	mean(subset(datMat, L1_Imp == 1 & L2_Imp == -1 )$phen),
	mean(subset(datMat, L1_Imp == -1 & L2_Imp == -1 )$phen)
	)
	
	)
	
	
	# boxplot(phen~.*.,data.frame(datMat))

	
	
	
	FILENAME<-paste(BARDIR,"/EpistatsisPlot","_",phenOfInterest,"_",DEgene,"_",ASEgene,"_",ImprintingScores$Marker.ID[locus1],"_",ImprintingScores$Marker.ID[locus2],".tiff", sep="")
	
	
	PLOT<-ggplot(DFMEANS,aes(x=L1_Imp, y=Means))+
	facet_grid(rows = vars(L2_Imp) )+
	geom_bar(stat='identity')+
	scale_x_discrete()+
	ggtitle( paste(phenOfInterest," ~ ",DEgene," X ", ASEgene,"\n",ImprintingScores$Marker.ID[locus1]," : ",ImprintingScores$Marker.ID[locus2], sep="") )+
	xlab("DE Genotype")+
	ylab("Phenotype")
	
	ggsave(as.character(FILENAME), PLOT, units="in",width=3.25,height=3)
	
	### Run ANOVA model and return F-value's 
	# c( summary(aov(phen~.*.,data.frame(datMat)))[[1]][["Pr(>F)"]][1:3],DEgene,ASEgene,phenOfInterest,ImprintingScores[locus1,4],ImprintingScores[locus2,4] )
	
}




PlotgenePairEpiTest<-function(gDE,gImp,SNPS,ImprintingScores,AdditiveScores,Phenotypes,phenOfInterest, BARDIR){

	### Pulls all loci for DE Gene
	gDELoci<-match(unique(subset(SNPS,Gene==gDE)$LocusID), ImprintingScores$LocusID)
	
	### Pulls all loci for ASE Gene
	gImpLoci<-match(unique(subset(SNPS,Gene== gImp)$LocusID), ImprintingScores$LocusID)
	
	### Generates matrix of All Loci Pairs
	LociMat<-matrix(unlist(lapply(gDELoci, function(A) lapply(gImpLoci,function(B) return(c(A,B))) )),byrow=TRUE,ncol=2)
	
	### Runs Test on all pairs in matrix. Returns list of F-statistics
	unlist((apply(LociMat,1,function(X) PlotImpTest(X[1], X[2], ImprintingScores, AdditiveScores, Phenotypes, phenOfInterest, gDE, gImp, BARDIR))))

}





}










## Running Epistasis Tests
{


PHENOTYPES<-sF16Phenotypes
PHENOTYPES<-F16Phenotypes

IMPRINTINGSCORES<-NetImprintingScores
ADDITIVESCORES<-NetAdditiveScores
MEDIANS<-NetMedian

### Subset network object to DE~ASE~Phen sets where both the ASE and DE gene has associated markers
Net<-Net[(Net$ASE.Gene %in% as.character(unique(GeneSNPs$Gene)) & Net$Target.Gene %in% as.character(unique(GeneSNPs$Gene))),]

### This is going to iterate through each line of the Net Data frame. For each row (N) it passes the DE name, ASE name, and Phenotype of interest, as well as ALL of the marker, imprinting scores, additive scores and phenotypes.

## Once it runs all the rows, it breaks down the object structure into a vectors and reassembles it into a matrix with 8 columns filling the matrix by rows. This is then convereted to a dataframe.

EpistasisData<-data.frame(
	matrix(
		unlist(
			c(
				lapply(1:nrow(Net), function(N)
					genePairEpiTest(
						as.character(Net[N,1]),
						as.character(Net[N,2]),
						GeneSNPs,
						IMPRINTINGSCORES,
						ADDITIVESCORES,
						PHENOTYPES,
						as.character(Net[N,5])
					)
				)
			)
		)
		,byrow=TRUE,ncol=8
	)
)




## Name columns of dataframe
colnames(EpistasisData)<-c("DE_Imp","ASE_Imp","DE_Imp:ASE_Imp","DE","ASE","Phenotype","DE Marker","ASE Marker")


EpistasisData$"DE_Imp:ASE_Imp"<-as.numeric(as.character(EpistasisData$"DE_Imp:ASE_Imp"))
EpistasisData$"DE_Imp"<-as.numeric(as.character(EpistasisData$"DE_Imp"))
EpistasisData$"ASE_Imp"<-as.numeric(as.character(EpistasisData$"ASE_Imp"))



}









### Defines a new object to be permuted
PERMSCORES<-IMPRINTINGSCORES

### Set number of iterations for permuting and define objects to hold permutation data
iterations<-600
IterOverall<-NULL
IterPreviousQuantile<-rep(0,101)
deltaQuants<-rep(1,iterations)


### I used a For loop here, which I would NOT normally do as it is ISANELY inefficient (literally exponentially slower), but coding it as an apply function would be much harder to read and do

### This will loop through all iterations

for(iter in 1:iterations){
	
	# ### For each marker it shuffles all of the genotypes WITH replacement.
	# for(Locus in 8:(ncol(PERMSCORES)-1)){
		# PERMSCORES[,Locus]<-sample(PERMSCORES[,Locus],nrow(PERMSCORES),replace=FALSE)
	# }
	
	
	### This is going to iterate through each line of the Net Data frame. For each row (N) it passes the DE name, ASE name, and Phenotype of interest, as well as ALL of the marker, PERMUTED imprinting scores, additive scores and phenotypes.

	## Once it runs all the rows, it breaks down the object structure into a vectors and reassembles it into a matrix with 8 columns filling the matrix by rows. This is then convereted to a dataframe.

	
	IterPermData<-data.frame(
		matrix(
			unlist(
				c(
					lapply(1:nrow(Net), function(N)
						genePairEpiTest(
						as.character(Net[N,1]),
						as.character(Net[N,2]),
						GeneSNPs,
						IMPRINTINGSCORES,
						ADDITIVESCORES,
						PHENOTYPES,
						as.character(Net[N,5]),
						TRUE
						)
					)
				)
			)
			,byrow=TRUE,ncol=8
		)
	)
	
	
	colnames(IterPermData)<-c("DE_Imp","ASE_Imp","DE_Imp:ASE_Imp","DE","ASE","Phenotype","DE Marker","ASE Marker")
	
	### Here it measures the 1-100th quantiles before the latest permutation data is added
	IterPreviousQuantile<-quantile(as.numeric(as.character(IterOverall$"DE_Imp:ASE_Imp")),seq(0,1,0.01), na.rm=TRUE)
	
	### Adding the newest permutation data
	IterOverall<-rbind(IterPermData,IterOverall)
	
	### Clean permutation data object to ensure data generated in the next cycle is not a repeat of this one
	IterPermData<-NULL
	
	### Measure the 1-100th quantiles after the latest permutation data is added
	IterCurrentQuantile<-quantile(as.numeric(as.character(IterOverall$"DE_Imp:ASE_Imp")),seq(0,1,0.01), na.rm=TRUE)
	
	### Calculate sum of absolute value of deviation of quantile values between this iteration and overally permuted data
	deltaQuants[iter]<-sum(abs(IterCurrentQuantile-IterPreviousQuantile))
	
}





### Save permutation data figure


Stabilizationdata<-data.frame(
DeltaQuants=deltaQuants[!is.na(deltaQuants)],
Iterations=2:length(deltaQuants)
)

StabilityPlot<-ggplot(Stabilizationdata,aes(x= Iterations, y= DeltaQuants))+
geom_point(size=0.1)+
ylab("Change in quantiles")

ggsave("NullDistributionStabilization.tiff",StabilityPlot,units="in",height=3,width=3)



IterOverall$"DE_Imp:ASE_Imp"<-as.numeric(as.character(IterOverall $"DE_Imp:ASE_Imp"))
IterOverall$"DE_Imp"<-as.numeric(as.character(IterOverall $"DE_Imp"))
IterOverall$"ASE_Imp"<-as.numeric(as.character(IterOverall $"ASE_Imp"))


write.table(IterOverall,file="PermutedData.txt",row.names=FALSE,quote=FALSE,sep="\t")


### Using the LSR package build empirical cumulative null distribution of F-statistics with ecdf function for each term
library(lsr)






EpistasisData<-EpistasisData[order(EpistasisData$Phenotype),]




FObj<-NULL

DEObj<-NULL
ASEObj<-NULL

qObj<-NULL
AppendedFObj<-NULL
AppendedDEObj<-NULL
AppendedASEObj<-NULL


for(PHEN in unique(EpistasisData$Phenotype)){
	
	print(PHEN)
	
	RealDist_DE.ASE<-ecdf( subset(EpistasisData, Phenotype==PHEN)$"DE_Imp:ASE_Imp" )
	NullDist_DE.ASE<-ecdf( subset(IterOverall, Phenotype==PHEN)$"DE_Imp:ASE_Imp" )
	
	RealDist_DE<-ecdf( subset(EpistasisData, Phenotype==PHEN)$"DE_Imp" )
	NullDist_DE<-ecdf( subset(IterOverall, Phenotype==PHEN)$"DE_Imp" )
	
	RealDist_ASE<-ecdf( subset(EpistasisData, Phenotype==PHEN)$"ASE_Imp" )
	NullDist_ASE<-ecdf( subset(IterOverall, Phenotype==PHEN)$"ASE_Imp" )
	
	

	
	
	qObj<-try(qvalue( subset(EpistasisData, Phenotype==PHEN)$"DE_Imp:ASE_Imp" ), TRUE)
	
	#print(qObj)
	
	
	
	
	

	FObj<-data.frame(
	NULLDat=unlist(lapply(subset(EpistasisData, Phenotype==PHEN)$"DE_Imp:ASE_Imp",function(t) NullDist_DE.ASE(t) )),
	REALDat=unlist(lapply(subset(EpistasisData, Phenotype==PHEN)$"DE_Imp:ASE_Imp",function(t) RealDist_DE.ASE(t) ))
	)



	DEObj<-data.frame(
	NULLDat=unlist(lapply(subset(EpistasisData, Phenotype==PHEN)$"DE_Imp",function(t) NullDist_DE(t) )),
	REALDat=unlist(lapply(subset(EpistasisData, Phenotype==PHEN)$"DE_Imp",function(t) RealDist_DE(t) ))
	)
	
	
	ASEObj<-data.frame(
	NULLDat=unlist(lapply(subset(EpistasisData, Phenotype==PHEN)$"ASE_Imp",function(t) NullDist_ASE(t) )),
	REALDat=unlist(lapply(subset(EpistasisData, Phenotype==PHEN)$"ASE_Imp",function(t) RealDist_ASE(t) ))
	)	
	
	
	
	FObj$Ratio<-(FObj$NULLDat+1e-16)/(FObj$REALDat+1e-16)
	
	DEObj$Ratio<-(DEObj$NULLDat+1e-16)/(DEObj$REALDat+1e-16)
	
	ASEObj$Ratio<-(ASEObj$NULLDat+1e-16)/(ASEObj$REALDat+1e-16)
	
	###	
	
	EPIST_MTC_PLOT<-plotPermFDRQC(subset(EpistasisData, Phenotype==PHEN)$"DE_Imp:ASE_Imp",FObj$Ratio,0.1,subset(IterOverall, Phenotype==PHEN)$"DE_Imp:ASE_Imp")
	ggsave(paste("Multiple_Tests_Correction_",PHEN,"_","ASExDE.tiff",sep=""), EPIST_MTC_PLOT)

	ASE_MTC_PLOT<-plotPermFDRQC(subset(EpistasisData, Phenotype==PHEN)$"ASE_Imp", ASEObj$Ratio,0.1,subset(IterOverall, Phenotype==PHEN)$"ASE_Imp")
	ggsave(paste("Multiple_Tests_Correction_",PHEN,"_","ASE.tiff",sep=""), ASE_MTC_PLOT)

	DE_MTC_PLOT<-plotPermFDRQC(subset(EpistasisData, Phenotype==PHEN)$"DE_Imp", DEObj$Ratio,0.1,subset(IterOverall, Phenotype==PHEN)$"DE_Imp")
	ggsave(paste("Multiple_Tests_Correction_",PHEN,"_","DE.tiff",sep=""), DE_MTC_PLOT)
	
	
	
	AppendedDEObj <-c(AppendedDEObj,DEObj$Ratio)
	AppendedASEObj <-c(AppendedASEObj,ASEObj$Ratio)
		
	
	if( is.na(sum(FObj$Ratio<=0.1)) ){
		print(FObj)
	} else {
		
		if(isFALSE(class(qObj)=="try-error")){
				print( c( sum(FObj$Ratio<=0.1), sum(qObj$qvalues<=0.1)) )
				AppendedFObj<-c(AppendedFObj,FObj$Ratio)
			} else {
				print( c( sum(FObj$Ratio<=0.1), NA) )
				AppendedFObj<-c(AppendedFObj,FObj$Ratio)
			}

	}
		
}



EpistasisData$DE.ASE_fdr<-AppendedFObj
EpistasisData$DE_fdr<-AppendedDEObj
EpistasisData$ASE_fdr<-AppendedASEObj


subset(EpistasisData, DE.ASE_fdr<=0.1)
subset(EpistasisData, DE_fdr<=0.1)
subset(EpistasisData, ASE_fdr<=0.1)




PlotgenePairEpiTest(
						"Car3",
						"Cdkn1c",
						GeneSNPs,
						IMPRINTINGSCORES,
						ADDITIVESCORES,
						PHENOTYPES,
						as.character("AUC2"),
						Bar_Dir
)





SIG_EPIST_NET<-subset(EpistasisData, DE.ASE_fdr<=0.1)[,c(4,5,3,6,6)]


lapply(1:nrow(SIG_EPIST_NET), function(N)
	PlotgenePairEpiTest(
		as.character(SIG_EPIST_NET[N,1]),
		as.character(SIG_EPIST_NET[N,2]),
		GeneSNPs,
		IMPRINTINGSCORES,
		ADDITIVESCORES,
		PHENOTYPES,
		as.character(SIG_EPIST_NET[N,5]),
		Bar_Dir
		)
				)






length(AppendedFObj)
dim(EpistasisData)




###### Only look at the closest marker


ClosestMarkersData<-NULL
ClosestAppendedFObj<-NULL

for(LINE in 1:nrow(Net)){
	
	GENE.ASE<-as.character(Net[LINE,2])
	GENE.DE<-as.character(Net[LINE,1])
	PHENO<-as.character(Net[LINE,5])
	
	ASE.Locus<-as.numeric(as.character(subset(NetMedian,Gene==GENE.ASE)$Locus))
	DE.Locus<-as.numeric(as.character(subset(NetMedian,Gene==GENE.DE)$Locus))
	
	
	PieceOfData<-subset(EpistasisData, DE == GENE.DE & ASE == GENE.ASE & Phenotype == PHENO)
	

	DEMARKERS<-unlist(lapply(as.character(unique(PieceOfData$"DE Marker")),function(DEMARKER) DEMARKER))

	DEDISTANCES<-unlist(lapply(as.character(unique(PieceOfData$"DE Marker")),function(DEMARKER) abs(DE.Locus-subset(NetImprintingScores, Marker.ID==DEMARKER)$Position.mm9.build37.)))

	Closest.DE<-DEMARKERS[DEDISTANCES==min(DEDISTANCES)]


	ASEMARKERS<-unlist(lapply(as.character(unique(PieceOfData$"ASE Marker")),function(ASEMARKER) ASEMARKER))

	ASEDISTANCES<-unlist(lapply(as.character(unique(PieceOfData$"ASE Marker")),function(ASEMARKER) abs(ASE.Locus-subset(NetImprintingScores, Marker.ID==ASEMARKER)$Position.mm9.build37.)))
	
	Closest.ASE<-ASEMARKERS[ASEDISTANCES==min(ASEDISTANCES)]
	
	
	ClosestMarkersData<-rbind(PieceOfData[PieceOfData$"DE Marker"== Closest.DE & PieceOfData$"ASE Marker"== Closest.ASE,],ClosestMarkersData)
		
}




ClosestMarkerIterOverall<-NULL

for(LINE in 1:nrow(ClosestMarkersData) ){	
	
	PossibleLines<-subset(IterOverall, DE==ClosestMarkersData[LINE,]$DE & ASE==ClosestMarkersData[LINE,]$ASE & Phenotype==ClosestMarkersData[LINE,]$Phenotype)
	
	ClosestMarkerIterOverall <- rbind( PossibleLines[PossibleLines$"DE Marker" == ClosestMarkersData[LINE,]$"DE Marker" & PossibleLines$"ASE Marker" == ClosestMarkersData[LINE,]$"ASE Marker",] , ClosestMarkerIterOverall )
	
}






FObj<-NULL
qObj<-NULL
ClosestAppendedDEObj<-NULL
ClosestAppendedASEObj<-NULL
ClosestAppendedFObj<-NULL

ClosestMarkersData <-ClosestMarkersData[order(ClosestMarkersData$Phenotype),]



for(PHEN in unique(ClosestMarkersData$Phenotype)){
	
	print(PHEN)
	
	RealDist_DE.ASE<-ecdf( subset(ClosestMarkersData, Phenotype==PHEN)$"DE_Imp:ASE_Imp" )
	NullDist_DE.ASE<-ecdf( subset(ClosestMarkerIterOverall, Phenotype==PHEN)$"DE_Imp:ASE_Imp" )
	
	
	RealDist_DE<-ecdf( subset(ClosestMarkersData, Phenotype==PHEN)$"DE_Imp" )
	NullDist_DE<-ecdf( subset(ClosestMarkerIterOverall, Phenotype==PHEN)$"DE_Imp" )
	
	RealDist_ASE<-ecdf( subset(ClosestMarkersData, Phenotype==PHEN)$"ASE_Imp" )
	NullDist_ASE<-ecdf( subset(ClosestMarkerIterOverall, Phenotype==PHEN)$"ASE_Imp" )

	
	
	qObj<-try(qvalue( subset(ClosestMarkersData, Phenotype==PHEN)$"DE_Imp:ASE_Imp" ), TRUE)
	
	#print(qObj)



	FObj<-data.frame(
	NULLDat=unlist(lapply(subset(ClosestMarkersData, Phenotype==PHEN)$"DE_Imp:ASE_Imp",function(t) NullDist_DE.ASE(t) )),
	REALDat=unlist(lapply(subset(ClosestMarkersData, Phenotype==PHEN)$"DE_Imp:ASE_Imp",function(t) RealDist_DE.ASE(t) ))
	)
	
	DEObj<-data.frame(
	NULLDat=unlist(lapply(subset(ClosestMarkersData, Phenotype==PHEN)$"DE_Imp",function(t) NullDist_DE(t) )),
	REALDat=unlist(lapply(subset(ClosestMarkersData, Phenotype==PHEN)$"DE_Imp",function(t) RealDist_DE(t) ))
	)
	
	ASEObj<-data.frame(
	NULLDat=unlist(lapply(subset(ClosestMarkersData, Phenotype==PHEN)$"ASE_Imp",function(t) NullDist_ASE(t) )),
	REALDat=unlist(lapply(subset(ClosestMarkersData, Phenotype==PHEN)$"ASE_Imp",function(t) RealDist_ASE(t) ))
	)	

	
	
	
	FObj$Ratio<-(FObj$NULLDat+1e-16)/(FObj$REALDat+1e-16)
	
	DEObj$Ratio<-(DEObj$NULLDat+1e-16)/(DEObj$REALDat+1e-16)
	
	ASEObj$Ratio<-(ASEObj$NULLDat+1e-16)/(ASEObj$REALDat+1e-16)
	
	
	EPIST_MTC_PLOT<-plotPermFDRQC(subset(ClosestMarkersData, Phenotype==PHEN)$"DE_Imp:ASE_Imp",FObj$Ratio,0.1,subset(ClosestMarkerIterOverall, Phenotype==PHEN)$"DE_Imp:ASE_Imp")
	ggsave(paste("ClosestMarker_Multiple_Tests_Correction_",PHEN,"_","ASExDE.tiff",sep=""), EPIST_MTC_PLOT)

	ASE_MTC_PLOT<-plotPermFDRQC(subset(ClosestMarkersData, Phenotype==PHEN)$"ASE_Imp", ASEObj$Ratio,0.1,subset(ClosestMarkerIterOverall, Phenotype==PHEN)$"ASE_Imp")
	ggsave(paste("ClosestMarker_Multiple_Tests_Correction_",PHEN,"_","ASE.tiff",sep=""), ASE_MTC_PLOT)

	DE_MTC_PLOT<-plotPermFDRQC(subset(ClosestMarkersData, Phenotype==PHEN)$"DE_Imp", DEObj$Ratio,0.1,subset(ClosestMarkerIterOverall, Phenotype==PHEN)$"DE_Imp")
	ggsave(paste("ClosestMarker_Multiple_Tests_Correction_",PHEN,"_","DE.tiff",sep=""), DE_MTC_PLOT)

	
	
	ClosestAppendedDEObj <-c(ClosestAppendedDEObj,DEObj$Ratio)
	ClosestAppendedASEObj <-c(ClosestAppendedASEObj,ASEObj$Ratio)
	
	
	if( is.na(sum(FObj$Ratio<=0.1)) ){
		print(FObj)
	} else {
		
		if(isFALSE(class(qObj)=="try-error")){
				print( c( sum(FObj$Ratio<=0.1), sum(qObj$qvalues<=0.1)) )
				ClosestAppendedFObj <-c(ClosestAppendedFObj,FObj$Ratio)

			} else {
				print( c( sum(FObj$Ratio<=0.1), NA) )
				ClosestAppendedFObj <-c(ClosestAppendedFObj,FObj$Ratio)

			}

	}
		
}

ClosestMarkersData$DE.ASE_fdr<-ClosestAppendedFObj

ClosestMarkersData$DE_fdr<-ClosestAppendedDEObj
ClosestMarkersData$ASE_fdr<-ClosestAppendedASEObj



subset(ClosestMarkersData, DE.ASE_fdr<=0.1)
subset(ClosestMarkersData, DE_fdr<=0.1)
subset(ClosestMarkersData, ASE_fdr<=0.1)


subset(EpistasisData, DE.ASE_fdr<=0.1)
subset(EpistasisData, DE_fdr<=0.1)
subset(EpistasisData, ASE_fdr<=0.1)








# # RealDist_DE<-ecdf(EpistasisData$DE_Imp)
# NullDist_DE<-ecdf(IterOverall$DE_Imp)

# RealDist_ASE<-ecdf(as.numeric(as.character(EpistasisData$"ASE_Imp")))
# NullDist_ASE<-ecdf(as.numeric(as.character(IterOverall$"ASE_Imp")))

# RealDist_DE.ASE<-ecdf(as.numeric(as.character(EpistasisData$"DE_Imp:ASE_Imp")))
# NullDist_DE.ASE<-ecdf(as.numeric(as.character(IterOverall$"DE_Imp:ASE_Imp")))




# # # unlist(lapply(EpistasisData$DE_Imp,function(t) NullDist_DE(t)))
# EpistasisData$DE_fdr<-unlist(lapply(EpistasisData$DE_Imp,function(t) NullDist_DE(t)/RealDist_DE(t) ))
# EpistasisData$ASE_fdr<-unlist(lapply(EpistasisData$ASE_Imp,function(t) NullDist_ASE(t)/RealDist_ASE(t) ))
# EpistasisData$DE.ASE_fdr<-unlist(lapply(EpistasisData$"DE_Imp:ASE_Imp",function(t) NullDist_DE.ASE(t)/RealDist_DE.ASE(t) ))



# ### For each real test, Calculate probabilty of observing an F-statistic at least as large under the respective NULL model 
# EpistasisData$DE_P<-unlist(lapply( as.numeric(as.character(EpistasisData$"DE_Imp")) ,function(t) 1-NullDist_DE(t)))
# EpistasisData$ASE_P<-unlist(lapply( as.numeric(as.character(EpistasisData$"ASE_Imp")),function(t) 1-NullDist_ASE(t)))
# EpistasisData$DE.ASE_P<-unlist(lapply( as.numeric(as.character(EpistasisData$"DE_Imp:ASE_Imp")),function(t) 1-NullDist_DE.ASE(t)))




# # hist(qvalue(EpistasisData$DE_Imp))

# hist(qvalue(EpistasisData$ASE_Imp))

# hist(qvalue(EpistasisData$"DE_Imp:ASE_Imp"))




# ### Return number of significant DE~ASE~Phenotype for each term
# sum(EpistasisData$DE_P<0.05,na.rm=TRUE)
# sum(EpistasisData$ASE_P<0.05,na.rm=TRUE)
# sum(EpistasisData$DE.ASE_P<0.05,na.rm=TRUE)


### Return number of significant DE~ASE~Phenotype for each term
sum(EpistasisData$DE_fdr<=0.1,na.rm=TRUE)
sum(EpistasisData$ASE_fdr<=0.1,na.rm=TRUE)
sum(EpistasisData$DE.ASE_fdr<=0.1,na.rm=TRUE)

sum(ClosestMarkersData$DE_fdr<=0.1,na.rm=TRUE)
sum(ClosestMarkersData$ASE_fdr<=0.1,na.rm=TRUE)
sum(ClosestMarkersData$DE.ASE_fdr<=0.1,na.rm=TRUE)


subset(ClosestMarkersData,DE.ASE_fdr<=0.1)
subset(EpistasisData,DE.ASE_fdr<=0.1)



### Return unqiue significant DE~ASE epistatic pairs
unique(EpistasisData[EpistasisData$DE.ASE_fdr<=0.1,4:5])

unique(ClosestMarkersData[ClosestMarkersData $DE.ASE_fdr<=0.1,4:5])



### Return number of unqiue DE and Number of unique ASE genes showing significant epistasis
length(unique(EpistasisData[EpistasisData$DE.ASE_fdr <=0.1,4]))
length(unique(EpistasisData[EpistasisData$DE.ASE_fdr <=0.1,5]))



### Save data
write.table(EpistasisData,file="AllEpistasisData.txt",row.names=FALSE,quote=FALSE,sep="\t")


### Save data
write.table(ClosestMarkersData,file="Closest_AllEpistasisData.txt",row.names=FALSE,quote=FALSE,sep="\t")





