set.seed(0)
rm(list=ls())

##Import Data
{

WATnetmm9<-read.delim("WATNet_mm9.BED",header=FALSE,sep="\t")
colnames(WATnetmm9)<-c("Chr","Start","Stop","ENSEMBL","Gene")

WATnetMedian<-data.frame(matrix(unlist( lapply(1:length(unique(WATnetmm9$Gene)), function(n) c( mean(c(min(unlist(subset(WATnetmm9,Gene==unique(WATnetmm9$Gene)[n])[2:3])),max(unlist(subset(WATnetmm9,Gene==unique(WATnetmm9$Gene)[n])[2:3])))), as.character(unique(WATnetmm9$Gene)[n]) ) )), byrow=TRUE, ncol=2))

colnames(WATnetMedian)<-c("Locus","Gene")

WATnetAdditiveScores<-read.csv("WATNetAdditiveScores.csv")
WATnetAdditiveScores$LocusID<-apply(WATnetAdditiveScores[,c(2,5)],1,paste,collapse="_")
WATnetAdditiveScores$Marker.ID<-as.character(WATnetAdditiveScores$Marker.ID)

WATnetImprintingScores<-read.csv("WATNetImprintingScores.csv")
WATnetImprintingScores$LocusID<-apply(WATnetImprintingScores[,c(2,5)],1,paste,collapse="_")
WATnetImprintingScores$Marker.ID<-as.character(WATnetImprintingScores$Marker.ID)

F16Phenotypes<-read.csv("F16Phenotypes.csv")

WATnet<-read.delim("WATNetTable.txt",header=TRUE,sep="\t")

GeneSNPs<-read.delim("F16_SNPS_within_3Mb_of_WATNet.BED",header=FALSE,sep="\t")
colnames(GeneSNPs)<-c("Chr","Start","Stop","ENSEMBL","Gene")

GeneSNPs$LocusID<-apply(GeneSNPs[1:2],1,paste,collapse="_")

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

sF16Phenotypes$GlucWt1[as.numeric(rownames(data.frame(residuals(aov(GlucWt1~Sex*Diet*+Litter_Size,sF16Phenotypes)))))]<-residuals(aov(GlucWt1~Sex*Diet*+Litter_Size,sF16Phenotypes))

sF16Phenotypes$GlucWt2[as.numeric(rownames(data.frame(residuals(aov(GlucWt2~Sex*Diet*+Litter_Size,sF16Phenotypes)))))]<-residuals(aov(GlucWt2~Sex*Diet*+Litter_Size,sF16Phenotypes))

sF16Phenotypes$T0_1[as.numeric(rownames(data.frame(residuals(aov(T0_1~Sex*Diet*+Litter_Size,sF16Phenotypes)))))]<-residuals(aov(T0_1~Sex*Diet*+Litter_Size,sF16Phenotypes))

sF16Phenotypes$T15_1[as.numeric(rownames(data.frame(residuals(aov(T15_1~Sex*Diet*+Litter_Size,sF16Phenotypes)))))]<-residuals(aov(T15_1~Sex*Diet*+Litter_Size,sF16Phenotypes))

sF16Phenotypes$T30_1[as.numeric(rownames(data.frame(residuals(aov(T30_1~Sex*Diet*+Litter_Size,sF16Phenotypes)))))]<-residuals(aov(T30_1~Sex*Diet*+Litter_Size,sF16Phenotypes))

sF16Phenotypes$T60_1[as.numeric(rownames(data.frame(residuals(aov(T60_1~Sex*Diet*+Litter_Size,sF16Phenotypes)))))]<-residuals(aov(T60_1~Sex*Diet*+Litter_Size,sF16Phenotypes))

sF16Phenotypes$T120_1[as.numeric(rownames(data.frame(residuals(aov(T120_1~Sex*Diet*+Litter_Size,sF16Phenotypes)))))]<-residuals(aov(T120_1~Sex*Diet*+Litter_Size,sF16Phenotypes))

sF16Phenotypes$AUC1[as.numeric(rownames(data.frame(residuals(aov(AUC1~Sex*Diet*+Litter_Size,sF16Phenotypes)))))]<-residuals(aov(AUC1~Sex*Diet*+Litter_Size,sF16Phenotypes))


sF16Phenotypes$T0_2[as.numeric(rownames(data.frame(residuals(aov(T0_2~Sex*Diet*+Litter_Size,sF16Phenotypes)))))]<-residuals(aov(T0_2~Sex*Diet*+Litter_Size,sF16Phenotypes))

sF16Phenotypes$T15_2[as.numeric(rownames(data.frame(residuals(aov(T15_2~Sex*Diet*+Litter_Size,sF16Phenotypes)))))]<-residuals(aov(T15_2~Sex*Diet*+Litter_Size,sF16Phenotypes))

sF16Phenotypes$T30_2[as.numeric(rownames(data.frame(residuals(aov(T30_2~Sex*Diet*+Litter_Size,sF16Phenotypes)))))]<-residuals(aov(T30_2~Sex*Diet*+Litter_Size,sF16Phenotypes))

sF16Phenotypes$T60_2[as.numeric(rownames(data.frame(residuals(aov(T60_2~Sex*Diet*+Litter_Size,sF16Phenotypes)))))]<-residuals(aov(T60_2~Sex*Diet*+Litter_Size,sF16Phenotypes))

sF16Phenotypes$T120_2[as.numeric(rownames(data.frame(residuals(aov(T120_2~Sex*Diet*+Litter_Size,sF16Phenotypes)))))]<-residuals(aov(T120_2~Sex*Diet*+Litter_Size,sF16Phenotypes))

sF16Phenotypes$AUC2[as.numeric(rownames(data.frame(residuals(aov(AUC2~Sex*Diet*+Litter_Size,sF16Phenotypes)))))]<-residuals(aov(AUC2~Sex*Diet*+Litter_Size,sF16Phenotypes))


sF16Phenotypes$NecWeight[as.numeric(rownames(data.frame(residuals(aov(NecWeight~Sex*Diet*+Litter_Size,sF16Phenotypes)))))]<-residuals(aov(NecWeight~Sex*Diet*+Litter_Size,sF16Phenotypes))

sF16Phenotypes$Heart[as.numeric(rownames(data.frame(residuals(aov(Heart~Sex*Diet*+Litter_Size,sF16Phenotypes)))))]<-residuals(aov(Heart~Sex*Diet*+Litter_Size,sF16Phenotypes))

sF16Phenotypes$L.Kid[as.numeric(rownames(data.frame(residuals(aov(L.Kid~Sex*Diet*+Litter_Size,sF16Phenotypes)))))]<-residuals(aov(L.Kid~Sex*Diet*+Litter_Size,sF16Phenotypes))

sF16Phenotypes$R.Kid[as.numeric(rownames(data.frame(residuals(aov(R.Kid~Sex*Diet*+Litter_Size,sF16Phenotypes)))))]<-residuals(aov(R.Kid~Sex*Diet*+Litter_Size,sF16Phenotypes))

sF16Phenotypes$Kidney[as.numeric(rownames(data.frame(residuals(aov(Kidney~Sex*Diet*+Litter_Size,sF16Phenotypes)))))]<-residuals(aov(Kidney~Sex*Diet*+Litter_Size,sF16Phenotypes))

sF16Phenotypes$Spleen[as.numeric(rownames(data.frame(residuals(aov(Spleen~Sex*Diet*+Litter_Size,sF16Phenotypes)))))]<-residuals(aov(Spleen~Sex*Diet*+Litter_Size,sF16Phenotypes))

sF16Phenotypes$Liver[as.numeric(rownames(data.frame(residuals(aov(Liver~Sex*Diet*+Litter_Size,sF16Phenotypes)))))]<-residuals(aov(Liver~Sex*Diet*+Litter_Size,sF16Phenotypes))

sF16Phenotypes$Go.dal[as.numeric(rownames(data.frame(residuals(aov(Go.dal~Sex*Diet*+Litter_Size,sF16Phenotypes)))))]<-residuals(aov(Go.dal~Sex*Diet*+Litter_Size,sF16Phenotypes))

sF16Phenotypes$Kid_Fat[as.numeric(rownames(data.frame(residuals(aov(Kid_Fat~Sex*Diet*+Litter_Size,sF16Phenotypes)))))]<-residuals(aov(Kid_Fat~Sex*Diet*+Litter_Size,sF16Phenotypes))

sF16Phenotypes$Mesen[as.numeric(rownames(data.frame(residuals(aov(Mesen~Sex*Diet*+Litter_Size,sF16Phenotypes)))))]<-residuals(aov(Mesen~Sex*Diet*+Litter_Size,sF16Phenotypes))

sF16Phenotypes$Ing[as.numeric(rownames(data.frame(residuals(aov(Ing~Sex*Diet*+Litter_Size,sF16Phenotypes)))))]<-residuals(aov(Ing~Sex*Diet*+Litter_Size,sF16Phenotypes))

sF16Phenotypes$FFA[as.numeric(rownames(data.frame(residuals(aov(FFA~Sex*Diet*+Litter_Size,sF16Phenotypes)))))]<-residuals(aov(FFA~Sex*Diet*+Litter_Size,sF16Phenotypes))

sF16Phenotypes$GLC[as.numeric(rownames(data.frame(residuals(aov(GLC~Sex*Diet*+Litter_Size,sF16Phenotypes)))))]<-residuals(aov(GLC~Sex*Diet*+Litter_Size,sF16Phenotypes))

sF16Phenotypes$CHOL[as.numeric(rownames(data.frame(residuals(aov(CHOL~Sex*Diet*+Litter_Size,sF16Phenotypes)))))]<-residuals(aov(CHOL~Sex*Diet*+Litter_Size,sF16Phenotypes))

sF16Phenotypes$TG[as.numeric(rownames(data.frame(residuals(aov(TG~Sex*Diet*+Litter_Size,sF16Phenotypes)))))]<-residuals(aov(TG~Sex*Diet*+Litter_Size,sF16Phenotypes))

sF16Phenotypes$Insulin[as.numeric(rownames(data.frame(residuals(aov(Insulin~Sex*Diet*+Litter_Size,sF16Phenotypes)))))]<-residuals(aov(Insulin~Sex*Diet*+Litter_Size,sF16Phenotypes))

sF16Phenotypes$Tot_Fat[as.numeric(rownames(data.frame(residuals(aov(Tot_Fat~Sex*Diet*+Litter_Size,sF16Phenotypes)))))]<-residuals(aov(Tot_Fat~Sex*Diet*+Litter_Size,sF16Phenotypes))

unique(apply(sF16Phenotypes[7:8],1,function(S) paste(S,collapse=".")))

cohortColors<-c("darkRed","darkBlue","black","red","blue")[as.numeric(as.factor(apply(sF16Phenotypes[7:8],1,function(S) paste(S,collapse="."))))]


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


### Epistasis Test Functions
{

AddAndImpTest<-function(locus1,locus2,ImprintingScores,AdditiveScores,Phenotypes,phenOfInterest,DEgene,ASEgene){
	matchedPhenotypes<-Phenotypes[match(colnames(ImprintingScores[locus1,8:ncol(ImprintingScores)]),Phenotypes$Animal),] #Should this be locus1?
	matchedPhenIndex<-match(phenOfInterest,colnames(matchedPhenotypes))
	
	datMat<-cbind(t(data.frame(ImprintingScores[c(locus1,locus2),8:ncol(ImprintingScores)])), t(data.frame(AdditiveScores[c(locus1,locus2),8:ncol(AdditiveScores)])), matchedPhenotypes[,matchedPhenIndex])


	datMat<-cbind(t(data.frame(ImprintingScores[c(locus1,locus2),8:ncol(ImprintingScores)])), matchedPhenotypes[,matchedPhenIndex])

	
	colnames(datMat)<-c("L1_Imp", "L2_Imp","phen")
	datMat<-data.frame(apply(datMat,2,as.numeric))
	
	datMat<-datMat[!is.na(datMat$phen),]

	summary(aov(phen~.*.,data.frame(datMat)))

	#c(summary(aov(phen~.*.,data.frame(datMat)))[[1]][["Pr(>F)"]][1:8],DEgene,ASEgene,phenOfInterest,ImprintingScores[locus1,4],ImprintingScores[locus2,4])

	c( summary(aov(phen~.*.,data.frame(datMat)))[[1]][["Pr(>F)"]][1:3],DEgene,ASEgene,phenOfInterest,ImprintingScores[locus1,4],ImprintingScores[locus2,4] )
}

genePairEpiTest<-function(gDE,gImp,SNPS,ImprintingScores,AdditiveScores,Phenotypes,phenOfInterest){
	gDELoci<-match(unique(subset(SNPS,Gene==gDE)$LocusID), ImprintingScores$LocusID)
	gImpLoci<-match(unique(subset(SNPS,Gene== gImp)$LocusID), ImprintingScores$LocusID)
	LociMat<-matrix(unlist(lapply(gDELoci, function(A) lapply(gImpLoci,function(B) return(c(A,B))) )),byrow=TRUE,ncol=2)
	unlist((apply(LociMat,1,function(X) AddAndImpTest(X[1], X[2], ImprintingScores, AdditiveScores, Phenotypes, phenOfInterest, gDE, gImp))))
}

}


### All Contexts
{

## Running Epistasis Tests
{

PHENOTYPES<-sF16Phenotypes
IMPRINTINGSCORES<-WATnetImprintingScores
ADDITIVESCORES<-WATnetAdditiveScores
MEDIANS<-WATnetMedian

WATnet<-WATnet[(WATnet$ASE.Gene %in% as.character(unique(GeneSNPs$Gene)) & WATnet$Target.Gene %in% as.character(unique(GeneSNPs$Gene))),]


#PairList<-subset(WATnet, ASE.Gene=="Nnat")[1:2,1:2]
PairList<-WATnet[1:2]


EpistasisData<-data.frame(
	matrix(
		unlist(
				lapply(1:length(colnames(PHENOTYPES)[c(10:11,13:41)]), function (P)
					lapply(1:nrow(PairList), function(X) unlist(
							c(
								genePairEpiTest(
									as.character(PairList[X,1]),
									as.character(PairList[X,2]),
									GeneSNPs,
									IMPRINTINGSCORES,
									ADDITIVESCORES,
									PHENOTYPES,
									colnames(PHENOTYPES)[c(10:11,13:41)][P]
								)
							)
						)
					)
				)
		),byrow=TRUE,ncol=8
	)
)


colnames(EpistasisData)<-c("DE_Imp","ASE_Imp","DE_Imp:ASE_Imp","DE","ASE","Phenotype","DE Marker","ASE Marker")

}

## Subsetting Closet Marker
{

EpistasisData$DElocus<-as.numeric(unlist(lapply(EpistasisData$DE, function(G) as.character(MEDIANS$Locus[as.character(MEDIANS$Gene)==G]))))
EpistasisData$ASElocus<-as.numeric(unlist(lapply(EpistasisData$ASE, function(G) as.character(MEDIANS$Locus[as.character(MEDIANS$Gene)==G]))))

EpistasisData$DEMarkerlocus<-as.numeric(unlist(lapply(EpistasisData$"DE Marker", function(L) as.character(IMPRINTINGSCORES$Position.mm9.build37.[as.character(IMPRINTINGSCORES$Marker.ID)==L])  )))

EpistasisData$ASEMarkerlocus<-as.numeric(unlist(lapply(EpistasisData$"ASE Marker", function(L) as.character(IMPRINTINGSCORES$Position.mm9.build37.[as.character(IMPRINTINGSCORES$Marker.ID)==L])  )))

EpistasisData$DEMarkerGeneDistance<-abs(EpistasisData$DElocus-EpistasisData$DEMarkerlocus)
EpistasisData$ASEMarkerGeneDistance<-abs(EpistasisData$ASElocus-EpistasisData$ASEMarkerlocus)


closestMarkers<-rbind(data.frame(matrix(unlist(lapply(unique(EpistasisData$DE), function(gene) as.character(unlist(subset(EpistasisData,DE==gene)[which.min(subset(EpistasisData,DE==gene)$DEMarkerGeneDistance),c(9,12)])))),byrow=TRUE,ncol=2))
,data.frame(matrix(unlist(lapply(unique(EpistasisData$ASE), function(gene) as.character(unlist(subset(EpistasisData,ASE==gene)[which.min(subset(EpistasisData,ASE==gene)$ASEMarkerGeneDistance),c(10,13)])))),byrow=TRUE,ncol=2))
)

colnames(closestMarkers)<-c("Gene","Marker")


closestMarkOnlyEpistasisData<-EpistasisData[
as.character(EpistasisData$"DE Marker")==as.character(closestMarkers$Marker[match(EpistasisData$DE,closestMarkers$Gene)]) & as.character(EpistasisData$"ASE Marker")==as.character(closestMarkers$Marker[match(EpistasisData$ASE,closestMarkers$Gene)])
,]

}

## Muliple Test Correction

{
EpistasisData$DE_Imp_FDR<-p.adjust(as.numeric(as.character(EpistasisData$DE_Imp)),method="fdr")
EpistasisData$ASE_Imp_FDR<-p.adjust(as.numeric(as.character(EpistasisData$ASE_Imp)),method="fdr")

EpistasisData$"DE_Imp:ASE_Imp_FDR"<-p.adjust(as.numeric(as.character(EpistasisData$"DE_Imp:ASE_Imp")),method="fdr")


closestMarkOnlyEpistasisData $DE_Imp_FDR<-p.adjust(as.numeric(as.character(closestMarkOnlyEpistasisData $DE_Imp)),method="fdr")
closestMarkOnlyEpistasisData $ASE_Imp_FDR<-p.adjust(as.numeric(as.character(closestMarkOnlyEpistasisData $ASE_Imp)),method="fdr")

closestMarkOnlyEpistasisData $"DE_Imp:ASE_Imp_FDR"<-p.adjust(as.numeric(as.character(closestMarkOnlyEpistasisData $"DE_Imp:ASE_Imp")),method="fdr")
}

}


write.table(EpistasisData,file="AllEpistasisData.txt",row.names=FALSE,quote=FALSE,sep="\t")
write.table(closestMarkOnlyEpistasisData,file="closestMarkOnlyEpistasisData.txt",row.names=FALSE,quote=FALSE,sep="\t")





### HF
{

## Running Epistasis Tests
{

PHENOTYPES<-subset(tF16Phenotypes, Sex=="F" & Diet=="H")[]
IMPRINTINGSCORES<-WATnetImprintingScores
ADDITIVESCORES<-WATnetAdditiveScores
MEDIANS<-WATnetMedian

WATnet<-WATnet[(WATnet$ASE.Gene %in% as.character(unique(GeneSNPs$Gene)) & WATnet$Target.Gene %in% as.character(unique(GeneSNPs$Gene))),]


#PairList<-subset(WATnet, ASE.Gene=="Nnat")[1:2,1:2]
PairList<-WATnet[1:2]


EpistasisData<-data.frame(
	matrix(
		unlist(
				lapply(1:length(colnames(PHENOTYPES)[c(10:11,13:41)]), function (P)
					lapply(1:nrow(PairList), function(X) unlist(
							c(
								genePairEpiTest(
									as.character(PairList[X,1]),
									as.character(PairList[X,2]),
									GeneSNPs,
									IMPRINTINGSCORES,
									ADDITIVESCORES,
									PHENOTYPES,
									colnames(PHENOTYPES)[c(10:11,13:41)][P]
								)
							)
						)
					)
				)
		),byrow=TRUE,ncol=8
	)
)


colnames(EpistasisData)<-c("DE_Imp","ASE_Imp","DE_Imp:ASE_Imp","DE","ASE","Phenotype","DE Marker","ASE Marker")

}

## Subsetting Closet Marker
{

EpistasisData$DElocus<-as.numeric(unlist(lapply(EpistasisData$DE, function(G) as.character(MEDIANS$Locus[as.character(MEDIANS$Gene)==G]))))
EpistasisData$ASElocus<-as.numeric(unlist(lapply(EpistasisData$ASE, function(G) as.character(MEDIANS$Locus[as.character(MEDIANS$Gene)==G]))))

EpistasisData$DEMarkerlocus<-as.numeric(unlist(lapply(EpistasisData$"DE Marker", function(L) as.character(IMPRINTINGSCORES$Position.mm9.build37.[as.character(IMPRINTINGSCORES$Marker.ID)==L])  )))

EpistasisData$ASEMarkerlocus<-as.numeric(unlist(lapply(EpistasisData$"ASE Marker", function(L) as.character(IMPRINTINGSCORES$Position.mm9.build37.[as.character(IMPRINTINGSCORES$Marker.ID)==L])  )))

EpistasisData$DEMarkerGeneDistance<-abs(EpistasisData$DElocus-EpistasisData$DEMarkerlocus)
EpistasisData$ASEMarkerGeneDistance<-abs(EpistasisData$ASElocus-EpistasisData$ASEMarkerlocus)


closestMarkers<-rbind(data.frame(matrix(unlist(lapply(unique(EpistasisData$DE), function(gene) as.character(unlist(subset(EpistasisData,DE==gene)[which.min(subset(EpistasisData,DE==gene)$DEMarkerGeneDistance),c(9,12)])))),byrow=TRUE,ncol=2))
,data.frame(matrix(unlist(lapply(unique(EpistasisData$ASE), function(gene) as.character(unlist(subset(EpistasisData,ASE==gene)[which.min(subset(EpistasisData,ASE==gene)$ASEMarkerGeneDistance),c(10,13)])))),byrow=TRUE,ncol=2))
)

colnames(closestMarkers)<-c("Gene","Marker")


closestMarkOnlyEpistasisData<-EpistasisData[
as.character(EpistasisData$"DE Marker")==as.character(closestMarkers$Marker[match(EpistasisData$DE,closestMarkers$Gene)]) & as.character(EpistasisData$"ASE Marker")==as.character(closestMarkers$Marker[match(EpistasisData$ASE,closestMarkers$Gene)])
,]

}

## Muliple Test Correction

{
EpistasisData$DE_Imp_FDR<-p.adjust(as.numeric(as.character(EpistasisData$DE_Imp)),method="fdr")
EpistasisData$ASE_Imp_FDR<-p.adjust(as.numeric(as.character(EpistasisData$ASE_Imp)),method="fdr")


EpistasisData$"DE_Imp:ASE_Imp_FDR"<-p.adjust(as.numeric(as.character(EpistasisData$"DE_Imp:ASE_Imp")),method="fdr")



closestMarkOnlyEpistasisData $DE_Imp_FDR<-p.adjust(as.numeric(as.character(closestMarkOnlyEpistasisData $DE_Imp)),method="fdr")
closestMarkOnlyEpistasisData $ASE_Imp_FDR<-p.adjust(as.numeric(as.character(closestMarkOnlyEpistasisData $ASE_Imp)),method="fdr")


closestMarkOnlyEpistasisData $"DE_Imp:ASE_Imp_FDR"<-p.adjust(as.numeric(as.character(closestMarkOnlyEpistasisData $"DE_Imp:ASE_Imp")),method="fdr")
}

}

write.table(EpistasisData,file="AllEpistasisDataHF.txt",row.names=FALSE,quote=FALSE,sep="\t")
write.table(closestMarkOnlyEpistasisData,file="closestMarkOnlyEpistasisDataHF.txt",row.names=FALSE,quote=FALSE,sep="\t")


### LF
{

## Running Epistasis Tests
{

PHENOTYPES<-subset(tF16Phenotypes, Sex=="F" & Diet=="L")[]
IMPRINTINGSCORES<-WATnetImprintingScores
ADDITIVESCORES<-WATnetAdditiveScores
MEDIANS<-WATnetMedian

WATnet<-WATnet[(WATnet$ASE.Gene %in% as.character(unique(GeneSNPs$Gene)) & WATnet$Target.Gene %in% as.character(unique(GeneSNPs$Gene))),]


#PairList<-subset(WATnet, ASE.Gene=="Nnat")[1:2,1:2]
PairList<-WATnet[1:2]


EpistasisData<-data.frame(
	matrix(
		unlist(
				lapply(1:length(colnames(PHENOTYPES)[c(10:11,13:41)]), function (P)
					lapply(1:nrow(PairList), function(X) unlist(
							c(
								genePairEpiTest(
									as.character(PairList[X,1]),
									as.character(PairList[X,2]),
									GeneSNPs,
									IMPRINTINGSCORES,
									ADDITIVESCORES,
									PHENOTYPES,
									colnames(PHENOTYPES)[c(10:11,13:41)][P]
								)
							)
						)
					)
				)
		),byrow=TRUE,ncol=8
	)
)


colnames(EpistasisData)<-c("DE_Imp","ASE_Imp","DE_Imp:ASE_Imp","DE","ASE","Phenotype","DE Marker","ASE Marker")

}

## Subsetting Closet Marker
{

EpistasisData$DElocus<-as.numeric(unlist(lapply(EpistasisData$DE, function(G) as.character(MEDIANS$Locus[as.character(MEDIANS$Gene)==G]))))
EpistasisData$ASElocus<-as.numeric(unlist(lapply(EpistasisData$ASE, function(G) as.character(MEDIANS$Locus[as.character(MEDIANS$Gene)==G]))))

EpistasisData$DEMarkerlocus<-as.numeric(unlist(lapply(EpistasisData$"DE Marker", function(L) as.character(IMPRINTINGSCORES$Position.mm9.build37.[as.character(IMPRINTINGSCORES$Marker.ID)==L])  )))

EpistasisData$ASEMarkerlocus<-as.numeric(unlist(lapply(EpistasisData$"ASE Marker", function(L) as.character(IMPRINTINGSCORES$Position.mm9.build37.[as.character(IMPRINTINGSCORES$Marker.ID)==L])  )))

EpistasisData$DEMarkerGeneDistance<-abs(EpistasisData$DElocus-EpistasisData$DEMarkerlocus)
EpistasisData$ASEMarkerGeneDistance<-abs(EpistasisData$ASElocus-EpistasisData$ASEMarkerlocus)


closestMarkers<-rbind(data.frame(matrix(unlist(lapply(unique(EpistasisData$DE), function(gene) as.character(unlist(subset(EpistasisData,DE==gene)[which.min(subset(EpistasisData,DE==gene)$DEMarkerGeneDistance),c(9,12)])))),byrow=TRUE,ncol=2))
,data.frame(matrix(unlist(lapply(unique(EpistasisData$ASE), function(gene) as.character(unlist(subset(EpistasisData,ASE==gene)[which.min(subset(EpistasisData,ASE==gene)$ASEMarkerGeneDistance),c(10,13)])))),byrow=TRUE,ncol=2))
)

colnames(closestMarkers)<-c("Gene","Marker")


closestMarkOnlyEpistasisData<-EpistasisData[
as.character(EpistasisData$"DE Marker")==as.character(closestMarkers$Marker[match(EpistasisData$DE,closestMarkers$Gene)]) & as.character(EpistasisData$"ASE Marker")==as.character(closestMarkers$Marker[match(EpistasisData$ASE,closestMarkers$Gene)])
,]

}

## Muliple Test Correction

{
EpistasisData$DE_Imp_FDR<-p.adjust(as.numeric(as.character(EpistasisData$DE_Imp)),method="fdr")
EpistasisData$ASE_Imp_FDR<-p.adjust(as.numeric(as.character(EpistasisData$ASE_Imp)),method="fdr")


EpistasisData$"DE_Imp:ASE_Imp_FDR"<-p.adjust(as.numeric(as.character(EpistasisData$"DE_Imp:ASE_Imp")),method="fdr")

closestMarkOnlyEpistasisData $DE_Imp_FDR<-p.adjust(as.numeric(as.character(closestMarkOnlyEpistasisData $DE_Imp)),method="fdr")
closestMarkOnlyEpistasisData $ASE_Imp_FDR<-p.adjust(as.numeric(as.character(closestMarkOnlyEpistasisData $ASE_Imp)),method="fdr")

closestMarkOnlyEpistasisData $"DE_Imp:ASE_Imp_FDR"<-p.adjust(as.numeric(as.character(closestMarkOnlyEpistasisData $"DE_Imp:ASE_Imp")),method="fdr")
}

}

write.table(EpistasisData,file="AllEpistasisDataLF.txt",row.names=FALSE,quote=FALSE,sep="\t")
write.table(closestMarkOnlyEpistasisData,file="closestMarkOnlyEpistasisDataLF.txt",row.names=FALSE,quote=FALSE,sep="\t")

### HM
{

## Running Epistasis Tests
{

PHENOTYPES<-subset(tF16Phenotypes, Sex=="M" & Diet=="H")[]
IMPRINTINGSCORES<-WATnetImprintingScores
ADDITIVESCORES<-WATnetAdditiveScores
MEDIANS<-WATnetMedian

WATnet<-WATnet[(WATnet$ASE.Gene %in% as.character(unique(GeneSNPs$Gene)) & WATnet$Target.Gene %in% as.character(unique(GeneSNPs$Gene))),]


#PairList<-subset(WATnet, ASE.Gene=="Nnat")[1:2,1:2]
PairList<-WATnet[1:2]


EpistasisData<-data.frame(
	matrix(
		unlist(
				lapply(1:length(colnames(PHENOTYPES)[c(10:11,13:41)]), function (P)
					lapply(1:nrow(PairList), function(X) unlist(
							c(
								genePairEpiTest(
									as.character(PairList[X,1]),
									as.character(PairList[X,2]),
									GeneSNPs,
									IMPRINTINGSCORES,
									ADDITIVESCORES,
									PHENOTYPES,
									colnames(PHENOTYPES)[c(10:11,13:41)][P]
								)
							)
						)
					)
				)
		),byrow=TRUE,ncol=8
	)
)


colnames(EpistasisData)<-c("DE_Imp","ASE_Imp","DE_Imp:ASE_Imp","DE","ASE","Phenotype","DE Marker","ASE Marker")

}

## Subsetting Closet Marker
{

EpistasisData$DElocus<-as.numeric(unlist(lapply(EpistasisData$DE, function(G) as.character(MEDIANS$Locus[as.character(MEDIANS$Gene)==G]))))
EpistasisData$ASElocus<-as.numeric(unlist(lapply(EpistasisData$ASE, function(G) as.character(MEDIANS$Locus[as.character(MEDIANS$Gene)==G]))))

EpistasisData$DEMarkerlocus<-as.numeric(unlist(lapply(EpistasisData$"DE Marker", function(L) as.character(IMPRINTINGSCORES$Position.mm9.build37.[as.character(IMPRINTINGSCORES$Marker.ID)==L])  )))

EpistasisData$ASEMarkerlocus<-as.numeric(unlist(lapply(EpistasisData$"ASE Marker", function(L) as.character(IMPRINTINGSCORES$Position.mm9.build37.[as.character(IMPRINTINGSCORES$Marker.ID)==L])  )))

EpistasisData$DEMarkerGeneDistance<-abs(EpistasisData$DElocus-EpistasisData$DEMarkerlocus)
EpistasisData$ASEMarkerGeneDistance<-abs(EpistasisData$ASElocus-EpistasisData$ASEMarkerlocus)


closestMarkers<-rbind(data.frame(matrix(unlist(lapply(unique(EpistasisData$DE), function(gene) as.character(unlist(subset(EpistasisData,DE==gene)[which.min(subset(EpistasisData,DE==gene)$DEMarkerGeneDistance),c(9,12)])))),byrow=TRUE,ncol=2))
,data.frame(matrix(unlist(lapply(unique(EpistasisData$ASE), function(gene) as.character(unlist(subset(EpistasisData,ASE==gene)[which.min(subset(EpistasisData,ASE==gene)$ASEMarkerGeneDistance),c(10,13)])))),byrow=TRUE,ncol=2))
)

colnames(closestMarkers)<-c("Gene","Marker")


closestMarkOnlyEpistasisData<-EpistasisData[
as.character(EpistasisData$"DE Marker")==as.character(closestMarkers$Marker[match(EpistasisData$DE,closestMarkers$Gene)]) & as.character(EpistasisData$"ASE Marker")==as.character(closestMarkers$Marker[match(EpistasisData$ASE,closestMarkers$Gene)])
,]

}

## Muliple Test Correction

{
EpistasisData$DE_Imp_FDR<-p.adjust(as.numeric(as.character(EpistasisData$DE_Imp)),method="fdr")
EpistasisData$ASE_Imp_FDR<-p.adjust(as.numeric(as.character(EpistasisData$ASE_Imp)),method="fdr")

EpistasisData$"DE_Imp:ASE_Imp_FDR"<-p.adjust(as.numeric(as.character(EpistasisData$"DE_Imp:ASE_Imp")),method="fdr")


closestMarkOnlyEpistasisData $DE_Imp_FDR<-p.adjust(as.numeric(as.character(closestMarkOnlyEpistasisData $DE_Imp)),method="fdr")
closestMarkOnlyEpistasisData $ASE_Imp_FDR<-p.adjust(as.numeric(as.character(closestMarkOnlyEpistasisData $ASE_Imp)),method="fdr")

closestMarkOnlyEpistasisData $"DE_Imp:ASE_Imp_FDR"<-p.adjust(as.numeric(as.character(closestMarkOnlyEpistasisData $"DE_Imp:ASE_Imp")),method="fdr")
}

}

write.table(EpistasisData,file="AllEpistasisDataHM.txt",row.names=FALSE,quote=FALSE,sep="\t")
write.table(closestMarkOnlyEpistasisData,file="closestMarkOnlyEpistasisDataHM.txt",row.names=FALSE,quote=FALSE,sep="\t")


### LF
{

## Running Epistasis Tests
{

PHENOTYPES<-subset(tF16Phenotypes, Sex=="M" & Diet=="L")[]
IMPRINTINGSCORES<-WATnetImprintingScores
ADDITIVESCORES<-WATnetAdditiveScores
MEDIANS<-WATnetMedian

WATnet<-WATnet[(WATnet$ASE.Gene %in% as.character(unique(GeneSNPs$Gene)) & WATnet$Target.Gene %in% as.character(unique(GeneSNPs$Gene))),]


#PairList<-subset(WATnet, ASE.Gene=="Nnat")[1:2,1:2]
PairList<-WATnet[1:2]


EpistasisData<-data.frame(
	matrix(
		unlist(
				lapply(1:length(colnames(PHENOTYPES)[c(10:11,13:41)]), function (P)
					lapply(1:nrow(PairList), function(X) unlist(
							c(
								genePairEpiTest(
									as.character(PairList[X,1]),
									as.character(PairList[X,2]),
									GeneSNPs,
									IMPRINTINGSCORES,
									ADDITIVESCORES,
									PHENOTYPES,
									colnames(PHENOTYPES)[c(10:11,13:41)][P]
								)
							)
						)
					)
				)
		),byrow=TRUE,ncol=8
	)
)


colnames(EpistasisData)<-c("DE_Imp","ASE_Imp","DE_Imp:ASE_Imp","DE","ASE","Phenotype","DE Marker","ASE Marker")

}

## Subsetting Closet Marker
{

EpistasisData$DElocus<-as.numeric(unlist(lapply(EpistasisData$DE, function(G) as.character(MEDIANS$Locus[as.character(MEDIANS$Gene)==G]))))
EpistasisData$ASElocus<-as.numeric(unlist(lapply(EpistasisData$ASE, function(G) as.character(MEDIANS$Locus[as.character(MEDIANS$Gene)==G]))))

EpistasisData$DEMarkerlocus<-as.numeric(unlist(lapply(EpistasisData$"DE Marker", function(L) as.character(IMPRINTINGSCORES$Position.mm9.build37.[as.character(IMPRINTINGSCORES$Marker.ID)==L])  )))

EpistasisData$ASEMarkerlocus<-as.numeric(unlist(lapply(EpistasisData$"ASE Marker", function(L) as.character(IMPRINTINGSCORES$Position.mm9.build37.[as.character(IMPRINTINGSCORES$Marker.ID)==L])  )))

EpistasisData$DEMarkerGeneDistance<-abs(EpistasisData$DElocus-EpistasisData$DEMarkerlocus)
EpistasisData$ASEMarkerGeneDistance<-abs(EpistasisData$ASElocus-EpistasisData$ASEMarkerlocus)


closestMarkers<-rbind(data.frame(matrix(unlist(lapply(unique(EpistasisData$DE), function(gene) as.character(unlist(subset(EpistasisData,DE==gene)[which.min(subset(EpistasisData,DE==gene)$DEMarkerGeneDistance),c(9,12)])))),byrow=TRUE,ncol=2))
,data.frame(matrix(unlist(lapply(unique(EpistasisData$ASE), function(gene) as.character(unlist(subset(EpistasisData,ASE==gene)[which.min(subset(EpistasisData,ASE==gene)$ASEMarkerGeneDistance),c(10,13)])))),byrow=TRUE,ncol=2))
)

colnames(closestMarkers)<-c("Gene","Marker")


closestMarkOnlyEpistasisData<-EpistasisData[
as.character(EpistasisData$"DE Marker")==as.character(closestMarkers$Marker[match(EpistasisData$DE,closestMarkers$Gene)]) & as.character(EpistasisData$"ASE Marker")==as.character(closestMarkers$Marker[match(EpistasisData$ASE,closestMarkers$Gene)])
,]

}

## Muliple Test Correction

{
EpistasisData$DE_Imp_FDR<-p.adjust(as.numeric(as.character(EpistasisData$DE_Imp)),method="fdr")
EpistasisData$ASE_Imp_FDR<-p.adjust(as.numeric(as.character(EpistasisData$ASE_Imp)),method="fdr")

EpistasisData$"DE_Imp:ASE_Imp_FDR"<-p.adjust(as.numeric(as.character(EpistasisData$"DE_Imp:ASE_Imp")),method="fdr")



closestMarkOnlyEpistasisData $DE_Imp_FDR<-p.adjust(as.numeric(as.character(closestMarkOnlyEpistasisData $DE_Imp)),method="fdr")
closestMarkOnlyEpistasisData $ASE_Imp_FDR<-p.adjust(as.numeric(as.character(closestMarkOnlyEpistasisData $ASE_Imp)),method="fdr")

closestMarkOnlyEpistasisData $"DE_Imp:ASE_Imp_FDR"<-p.adjust(as.numeric(as.character(closestMarkOnlyEpistasisData $"DE_Imp:ASE_Imp")),method="fdr")
}

}

write.table(EpistasisData,file="AllEpistasisDataLM.txt",row.names=FALSE,quote=FALSE,sep="\t")
write.table(closestMarkOnlyEpistasisData,file="closestMarkOnlyEpistasisDataLM.txt",row.names=FALSE,quote=FALSE,sep="\t")





