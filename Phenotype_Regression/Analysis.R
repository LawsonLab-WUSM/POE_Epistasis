

rm(list=ls())
set.seed(0)

library(boot)
library(gplots)
library(RColorBrewer)
library(viridis)
library(ggplot2)
library(cowplot)


SAMPLE_INFORMATION_FILE<-"Sample_Info.txt"
PHENOTYPES_FILE<-"Phenotypes.txt"
NETWORK_FILE<-"filtEpi_FDR.txt"
POEQTL_NETWORK_FILE<-"poefiltEpiData.tsv"
NORMALIZED_READ_COUNTS_FILE<-"WATCounts.txt"
NCBI_ID_CONVERSION_FILE<-"IDconversions.txt"
LGSM_ID_CONVERSION_FILE<-"IDconversions_LGSM.txt"

NETWORK_TABLE_FILE<-'WATNetTable.tsv'




plotTransform<-function(Phenotype,Original,Transformed){
	png(filename=paste(Phenotype,"_TransformationPlot.png",sep=""), width=500,height=500,type=c("cairo"))
	par(mfrow=c(2,2))
	
	hist(Original,ylab="Frequency",xlab=Phenotype,main="Original")
	hist(Transformed,ylab="Frequency",xlab=Phenotype,main="Transformed")
	
	qqnorm(Original,pch=16,cex=0.5,main="Original")
	qqline(Original,col="red")
	
	qqnorm(Transformed,pch=16,cex=0.5,main="Tranformed")
	qqline(Transformed,col="red")
	dev.off()
}


rawPhenotypes<-read.delim(PHENOTYPES_FILE,header=TRUE)


{


cohortColors<-c("darkRed","red","black","darkBlue","blue")[as.numeric(as.factor(apply(rawPhenotypes[2:3],1,function(S) paste(S,collapse="."))))]


##### Transforming Phenotypes!

colnames(rawPhenotypes)[10:ncol(rawPhenotypes)][unlist(apply(rawPhenotypes[10:ncol(rawPhenotypes)],2,function(p) shapiro.test(p)[2]))<0.05]

unlist(apply(rawPhenotypes[10:ncol(rawPhenotypes)],2,function(p) shapiro.test(p)[2]))[unlist(apply(rawPhenotypes[10:ncol(rawPhenotypes)],2,function(p) shapiro.test(p)[2]))<0.05]


UNTRANS<-rawPhenotypes$Necweight
summary(aov(Necweight~cross+sex+diet,rawPhenotypes))
hist(rawPhenotypes$Necweight)
rawPhenotypes$Necweight[rawPhenotypes$sex=="F" & rawPhenotypes$diet=="H"]<-scale(subset(rawPhenotypes, sex=="F" & diet=="H")$Necweight)
rawPhenotypes$Necweight[rawPhenotypes$sex=="F" & rawPhenotypes$diet=="L"]<-scale(subset(rawPhenotypes, sex=="F" & diet=="L")$Necweight)
rawPhenotypes$Necweight[rawPhenotypes$sex=="M" & rawPhenotypes$diet=="H"]<-scale(subset(rawPhenotypes, sex=="M" & diet=="H")$Necweight)
rawPhenotypes$Necweight[rawPhenotypes$sex=="M" & rawPhenotypes$diet=="L"]<-scale(subset(rawPhenotypes, sex=="M" & diet=="L")$Necweight)
hist(rawPhenotypes$Necweight)
shapiro.test(rawPhenotypes$Necweight)

plotTransform("Necweight", UNTRANS,rawPhenotypes$Necweight)

png(filename=paste("Necweight","_ResidualPlot.png",sep=""), width=500,height=500,type=c("cairo"))
plot(rawPhenotypes$Necweight, UNTRANS,pch=16,cex=1,col=cohortColors,ylab="Raw Data",xlab="Residualized Transformed Data")
dev.off()


UNTRANS<-rawPhenotypes$GTTweight
summary(aov(GTTweight~cross+sex+diet,rawPhenotypes))
hist(rawPhenotypes$GTTweight)
rawPhenotypes$GTTweight[rawPhenotypes$sex=="F" & rawPhenotypes$diet=="H"]<-scale(subset(rawPhenotypes, sex=="F" & diet=="H")$GTTweight)
rawPhenotypes$GTTweight[rawPhenotypes$sex=="F" & rawPhenotypes$diet=="L"]<-scale(subset(rawPhenotypes, sex=="F" & diet=="L")$GTTweight)
rawPhenotypes$GTTweight[rawPhenotypes$sex=="M" & rawPhenotypes$diet=="H"]<-scale(subset(rawPhenotypes, sex=="M" & diet=="H")$GTTweight)
rawPhenotypes$GTTweight[rawPhenotypes$sex=="M" & rawPhenotypes$diet=="L"]<-scale(subset(rawPhenotypes, sex=="M" & diet=="L")$GTTweight)
hist(rawPhenotypes$GTTweight)
shapiro.test(rawPhenotypes$GTTweight)
plotTransform("GTTweight", UNTRANS,rawPhenotypes$GTTweight)

png(filename=paste("GTTweight","_ResidualPlot.png",sep=""), width=500,height=500,type=c("cairo"))
plot(rawPhenotypes$GTTweight, UNTRANS,pch=16,cex=1,col=cohortColors,ylab="Raw Data",xlab="Residualized Transformed Data")
dev.off()



UNTRANS<-rawPhenotypes$Week19
summary(aov(Week19~cross+sex+diet,rawPhenotypes))
hist(rawPhenotypes$Week19)
rawPhenotypes$Week19[rawPhenotypes$sex=="F" & rawPhenotypes$diet=="H"]<-scale(subset(rawPhenotypes, sex=="F" & diet=="H")$Week19)
rawPhenotypes$Week19[rawPhenotypes$sex=="F" & rawPhenotypes$diet=="L"]<-scale(subset(rawPhenotypes, sex=="F" & diet=="L")$Week19)
rawPhenotypes$Week19[rawPhenotypes$sex=="M" & rawPhenotypes$diet=="H"]<-scale(subset(rawPhenotypes, sex=="M" & diet=="H")$Week19)
rawPhenotypes$Week19[rawPhenotypes$sex=="M" & rawPhenotypes$diet=="L"]<-scale(subset(rawPhenotypes, sex=="M" & diet=="L")$Week19)
hist(rawPhenotypes$Week19)
shapiro.test(rawPhenotypes$Week19)
plotTransform("Week19", UNTRANS,rawPhenotypes$Week19)

png(filename=paste("Week19","_ResidualPlot.png",sep=""), width=500,height=500,type=c("cairo"))
plot(rawPhenotypes$Week19, UNTRANS,pch=16,cex=1,col=cohortColors,ylab="Raw Data",xlab="Residualized Transformed Data")
dev.off()




UNTRANS<-rawPhenotypes$averageLean
summary(aov(averageLean~cross+sex+diet,rawPhenotypes))
hist(rawPhenotypes$averageLean)
rawPhenotypes$averageLean[rawPhenotypes$sex=="F"]<-scale(subset(rawPhenotypes, sex=="F")$averageLean)
rawPhenotypes$averageLean[rawPhenotypes$sex=="M"]<-scale(subset(rawPhenotypes, sex=="M")$averageLean)
hist(rawPhenotypes$averageLean)
shapiro.test(rawPhenotypes$averageLean)
plotTransform("averageLean", UNTRANS,rawPhenotypes$averageLean)

png(filename=paste("averageLean","_ResidualPlot.png",sep=""), width=500,height=500,type=c("cairo"))
plot(rawPhenotypes$averageLean, UNTRANS,pch=16,cex=1,col=cohortColors,ylab="Raw Data",xlab="Residualized Transformed Data")
dev.off()



UNTRANS<-rawPhenotypes$FatLeanRatio
summary(aov(FatLeanRatio~cross+sex+diet,rawPhenotypes))
hist(rawPhenotypes$FatLeanRatio)
rawPhenotypes$FatLeanRatio[rawPhenotypes$sex=="F" & rawPhenotypes$diet=="H"]<-scale(subset(rawPhenotypes, sex=="F" & diet=="H")$FatLeanRatio)
rawPhenotypes$FatLeanRatio[rawPhenotypes$sex=="F" & rawPhenotypes$diet=="L"]<-scale(subset(rawPhenotypes, sex=="F" & diet=="L")$FatLeanRatio)
rawPhenotypes$FatLeanRatio[rawPhenotypes$sex=="M" & rawPhenotypes$diet=="H"]<-scale(subset(rawPhenotypes, sex=="M" & diet=="H")$FatLeanRatio)
rawPhenotypes$FatLeanRatio[rawPhenotypes$sex=="M" & rawPhenotypes$diet=="L"]<-scale(subset(rawPhenotypes, sex=="M" & diet=="L")$FatLeanRatio)
hist(rawPhenotypes$FatLeanRatio)
shapiro.test(rawPhenotypes$FatLeanRatio)
plotTransform("FatLeanRatio", UNTRANS,rawPhenotypes$FatLeanRatio)


png(filename=paste("FatLeanRatio","_ResidualPlot.png",sep=""), width=500,height=500,type=c("cairo"))
plot(rawPhenotypes$FatLeanRatio, UNTRANS,pch=16,cex=1,col=cohortColors,ylab="Raw Data",xlab="Residualized Transformed Data")
dev.off()



UNTRANS<-rawPhenotypes$Insulin

hist(rawPhenotypes$Insulin)
hist(log2(rawPhenotypes$Insulin))
shapiro.test(log2(rawPhenotypes$Insulin))
rawPhenotypes$Insulin<-log2(rawPhenotypes$Insulin)

plotTransform("Insulin", UNTRANS,rawPhenotypes$Insulin)

UNTRANS<-rawPhenotypes$FFA
hist(rawPhenotypes$FFA)
hist(log10(rawPhenotypes$FFA))
shapiro.test((rawPhenotypes$FFA))
shapiro.test(log10(rawPhenotypes$FFA))
rawPhenotypes$FFA<-log10(rawPhenotypes$FFA)

plotTransform("FFA",UNTRANS,rawPhenotypes$FFA)


UNTRANS<-rawPhenotypes$Cholesterol
hist(rawPhenotypes$Cholesterol)
hist(log2(rawPhenotypes$Cholesterol))
shapiro.test(log10(rawPhenotypes$Cholesterol))
rawPhenotypes$Cholesterol<-log10(rawPhenotypes$Cholesterol)
plotTransform("Cholesterol",UNTRANS,rawPhenotypes$Cholesterol)


UNTRANS<-rawPhenotypes$AUC.G.
hist(rawPhenotypes$AUC.G.)
hist(log2(rawPhenotypes$AUC.G.))
shapiro.test(log2(rawPhenotypes$AUC.G.))
rawPhenotypes$AUC.G.<-log2(rawPhenotypes$AUC.G.)

plotTransform("AUC.G.",UNTRANS,rawPhenotypes$AUC.G.)


	
}



Info<-read.delim(SAMPLE_INFORMATION_FILE,sep="\t",header=TRUE)


EpiNet<-read.delim(NETWORK_FILE,header=TRUE)
POEEpiNet<-read.delim(POEQTL_NETWORK_FILE,header=TRUE)



dim(EpiNet)

length(unique(EpiNet$gA))
length(unique(EpiNet$gB))

length(unique(c(as.character(unique(EpiNet$gA)),as.character(unique(EpiNet$gB)))))



Counts<-read.delim(NORMALIZED_READ_COUNTS_FILE,header=TRUE)



BadIDs<-c("ENSMUSG00000053935","ENSMUSG00000073617","ENSMUSG00000051586","ENSMUSG00000090399","ENSMUSG00000095832")

IDCon<-read.delim(NCBI_ID_CONVERSION_FILE)
matchedIDs<-subset(IDCon,toupper(Gene.stable.ID) %in% toupper(Counts$Gene) & !is.na(NCBI.gene.ID))

# IDCon<-read.delim(LGSM_ID_CONVERSION_FILE)
# IDCon<-subset(IDCon, !(Gene.stable.ID %in% BadIDs))
# matchedIDs<-subset(IDCon,(Gene.stable.ID) %in%(Counts$Gene))


Counts<-subset(Counts, Gene %in% as.character(matchedIDs$Gene.stable.ID))

Counts$Gene<-as.character(matchedIDs[match(toupper(Counts$Gene), toupper(matchedIDs$Gene.stable.ID)),]$Gene.name)
rownames(Counts)<-Counts$Gene
Counts<-Counts[,2:ncol(Counts)]


BarcodePos<-match(colnames(Counts[,1:ncol(Counts)]),Info$GTAC.Tag)

Cross<-Info$Cross[BarcodePos]
Diet<-Info$Diet[BarcodePos]
Sex<-Info$Sex[BarcodePos]


#### Normalize Functions
{

checkGeneNormality<-function(GENE){
	as.numeric(c(unlist(shapiro.test(as.numeric(Counts[rownames(Counts)==GENE,])))[2],
		summary(aov(as.numeric(Counts[rownames(Counts)==GENE,])~Sex+Diet+Cross, ))[[1]][["Pr(>F)"]][1:2]
	))
}

scaleGeneSex<-function(GENE){
	Counts[rownames(Counts)==GENE,Sex=="F"]<-scale(as.numeric(Counts[rownames(Counts)==GENE,Sex=="F"]))[,1]
	Counts[rownames(Counts)==GENE,Sex=="M"]<-scale(as.numeric(Counts[rownames(Counts)==GENE,Sex=="M"]))[,1]
	return(Counts[rownames(Counts)==GENE,])
}

scaleGeneDiet<-function(GENE){
	Counts[rownames(Counts)==GENE,Diet=="H"]<-scale(as.numeric(Counts[rownames(Counts)==GENE,Diet=="H"]))[,1]
	Counts[rownames(Counts)==GENE,Diet=="L"]<-scale(as.numeric(Counts[rownames(Counts)==GENE,Diet=="L"]))[,1]
	return(Counts[rownames(Counts)==GENE,])
}

scaleGeneSexDiet<-function(GENE){
	Counts[rownames(Counts)==GENE,Sex=="F" & Diet=="H"]<-scale(as.numeric(Counts[rownames(Counts)==GENE,Sex=="F" & Diet=="H"]))[,1]
	Counts[rownames(Counts)==GENE,Sex=="M" & Diet=="H"]<-scale(as.numeric(Counts[rownames(Counts)==GENE,Sex=="M" & Diet=="H"]))[,1]
	Counts[rownames(Counts)==GENE,Sex=="F" & Diet=="L"]<-scale(as.numeric(Counts[rownames(Counts)==GENE,Sex=="F" & Diet=="L"]))[,1]
	Counts[rownames(Counts)==GENE,Sex=="M" & Diet=="L"]<-scale(as.numeric(Counts[rownames(Counts)==GENE,Sex=="M" & Diet=="L"]))[,1]
	return(Counts[rownames(Counts)==GENE,])
}

checkAndNormalize<-function(GENE){
	PreTransformTests<-checkGeneNormality(GENE)
	TransformedData<-NULL
	if(PreTransformTests[1]<=0.05 & PreTransformTests[2]<=0.1 & PreTransformTests[3] <=0.1) TransformedData<-scaleGeneSexDiet(GENE)
	if(PreTransformTests[1]<=0.05 & PreTransformTests[2]>0.1 & PreTransformTests[3] <=0.1) TransformedData<-scaleGeneDiet(GENE)
	if(PreTransformTests[1]<=0.05 & PreTransformTests[2]<=0.1 & PreTransformTests[3] >0.1) TransformedData<-scaleGeneSex(GENE)
	if(PreTransformTests[1]<=0.05 & PreTransformTests[2]>0.1 & PreTransformTests[3] >0.1) TransformedData<-Counts[rownames(Counts)==GENE,]
	if(PreTransformTests[1]>0.05) TransformedData<-Counts[rownames(Counts)==GENE,]
	return(TransformedData)
}

}



geneList<-unique(c(as.character(unique(EpiNet$gA)),as.character(unique(EpiNet$gB))))[unique(c(as.character(unique(EpiNet$gA)),as.character(unique(EpiNet$gB)))) %in% rownames(Counts)]





Counts<-Counts[rownames(Counts) %in% geneList,]

for (GENE in as.character(rownames(Counts))){
	Counts[rownames(Counts)==GENE,]<-checkAndNormalize(GENE)
}



POEQTLgeneList<-unique(c(as.character(unique(POEEpiNet$gA)),as.character(unique(POEEpiNet$gB))))


# PhenotypeColumns<-c(10,11,12,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29)
PhenotypeColumns<-c(10,11,12,15,19,20,26,29)


CorMat<-NULL
CorMat<-matrix(
		byrow=TRUE,
		ncol=length(PhenotypeColumns),
		unlist(
			lapply(
				geneList,
				function(G) cor(cbind(
					as.numeric(Counts[rownames(Counts)==G,]),
					rawPhenotypes[match(colnames(Counts),rawPhenotypes$GTAC.Tag), PhenotypeColumns]
					),use="na.or.complete")[1,2:(length(PhenotypeColumns)+1)]
			)
		)
)




CorTestMat<-matrix(
		byrow=TRUE,
		ncol=length(PhenotypeColumns),
		unlist(
			lapply(
				geneList, function(G) lapply(PhenotypeColumns, function (P) cor.test(
				y=as.numeric(Counts[rownames(Counts)==G,]),
				x=rawPhenotypes[match(colnames(Counts),rawPhenotypes$GTAC.Tag),P]
				)[3]
				)

			)
		)
)


dim(CorMat)
dim(CorTestMat)

colnames(CorMat)<-colnames(rawPhenotypes)[PhenotypeColumns]
rownames(CorMat)<-geneList
CorMat[is.na(CorMat)]<-0


colnames(CorTestMat)<-colnames(rawPhenotypes)[PhenotypeColumns]
rownames(CorTestMat)<-geneList



library(qvalue)



QObj<-qvalue(as.data.frame(as.table(CorTestMat))[,3])

plot(QObj)
hist(QObj)


png(filename="Multiple_Tests_Correction.png", width=500,height=500,type=c("cairo"))
plot(QObj)
dev.off()



png(filename="Pvalue_Distribution.png", width=500,height=500,type=c("cairo"))
hist(QObj)
dev.off()


CorTestFDR<-matrix(QObj$qvalues,byrow=FALSE,ncol=length(PhenotypeColumns))


colnames(CorTestFDR)<-colnames(CorTestMat)
rownames(CorTestFDR)<-rownames(CorTestMat)



plot(y=-log10(CorTestFDR),c(CorMat),pch=16,cex=0.5,xlab="Pearson's Correlation",ylab="-Log10(P-value)",xlim=c(-1,1))
lines(c(-1.1,1.1),c(-log10(0.05),-log10(0.05)),col="grey",lwd=2)
lines(c(-0.5,-0.5),c(-999,999),col="red",lwd=2)
lines(c(0.5,0.5),c(-999,999),col="red",lwd=2)





PearsonsCutoff<-0.5

FDRcutoff<-0.05





AllCorTable<-as.data.frame(as.table(CorTestFDR))
colnames(AllCorTable)<-c("Gene","Phenotype","FDR")
AllCorTable$Pearsons<-as.data.frame(as.table(CorMat))[,3]


AllVolcanoPlot<-ggplot(AllCorTable,aes(x=Pearsons,y=-log10(FDR+1e-16) ))+
geom_point()+
geom_hline(yintercept=-log10(FDRcutoff), color="red",size=1)+
geom_vline(xintercept=PearsonsCutoff, color="blue", size=1)+
geom_vline(xintercept=-PearsonsCutoff, color="blue", size=1)+
ylab("-log10( FDR )")+
xlab("Pearson's correlation")+
xlim(c(-1,1))

ggsave("AllVolcanoPlot.tiff", AllVolcanoPlot)






GenesCorToTraits<-rownames(CorTestFDR)[apply(abs(CorTestFDR),1,min)<=FDRcutoff]


mainCorMat<-CorMat[rownames(CorMat)%in%GenesCorToTraits,]
mainCorMatFDR<-CorTestFDR[rownames(CorMat)%in%GenesCorToTraits,]

POEQTLmainCorMat<-CorMat[rownames(CorMat)%in%GenesCorToTraits & rownames(CorMat)%in%POEQTLgeneList,]
POEQTLmainCorMatFDR<-CorTestFDR[rownames(CorMat)%in%GenesCorToTraits & rownames(CorMat)%in%POEQTLgeneList,]


mainCorMat.Phen.dendro<-as.dendrogram(hclust(d = dist(x = t(mainCorMat) )))
mainCorMat.Phen.order<-order.dendrogram(mainCorMat.Phen.dendro)

mainCorMat.Gene.dendro<-as.dendrogram(hclust(d = dist(x = (mainCorMat) )))
mainCorMat.Gene.order<-order.dendrogram(mainCorMat.Gene.dendro)

mainCorMat<-mainCorMat[mainCorMat.Gene.order, mainCorMat.Phen.order]
mainCorMatFDR <-mainCorMatFDR[mainCorMat.Gene.order, mainCorMat.Phen.order]




POEQTLmainCorMat.Phen.dendro<-as.dendrogram(hclust(d = dist(x = t(POEQTLmainCorMat) )))
POEQTLmainCorMat.Phen.order<-order.dendrogram(POEQTLmainCorMat.Phen.dendro)

POEQTLmainCorMat.Gene.dendro<-as.dendrogram(hclust(d = dist(x = (POEQTLmainCorMat) )))
POEQTLmainCorMat.Gene.order<-order.dendrogram(POEQTLmainCorMat.Gene.dendro)

POEQTLmainCorMat<-POEQTLmainCorMat[POEQTLmainCorMat.Gene.order, POEQTLmainCorMat.Phen.order]
POEQTLmainCorMatFDR <-POEQTLmainCorMatFDR[POEQTLmainCorMat.Gene.order, POEQTLmainCorMat.Phen.order]








write.table(mainCorMat, file='NetCorMat.tsv', quote=FALSE, sep='\t', row.names=TRUE)
write.table(mainCorMatFDR, file='NetCorMatFDR.tsv', quote=FALSE, sep='\t', row.names=TRUE)

write.table(POEQTLmainCorMat, file='POEQTLNetCorMat.tsv', quote=FALSE, sep='\t', row.names=TRUE)
write.table(POEQTLmainCorMatFDR, file='POEQTLNetCorMatFDR.tsv', quote=FALSE, sep='\t', row.names=TRUE)




mainCorTable<-as.data.frame(as.table(mainCorMatFDR))
colnames(mainCorTable)<-c("Gene","Phenotype","FDR")
mainCorTable$Pearsons<-as.data.frame(as.table(mainCorMat))[,3]




write.table(mainCorTable,file="NetCorTable.tsv", quote=FALSE, sep='\t', row.names=FALSE)



POEQTLmainCorTable<-as.data.frame(as.table(POEQTLmainCorMatFDR))
colnames(POEQTLmainCorTable)<-c("Gene","Phenotype","FDR")
POEQTLmainCorTable$Pearsons<-as.data.frame(as.table(POEQTLmainCorMat))[,3]


write.table(POEQTLmainCorTable,file="POEQTLNetTable.tsv", quote=FALSE, sep='\t', row.names=FALSE)

filtmainCorTable<-subset(mainCorTable, abs(Pearsons)>=PearsonsCutoff & FDR<=FDRcutoff)

POEQTLfiltmainCorTable<-subset(POEQTLmainCorTable, abs(Pearsons)>=PearsonsCutoff & FDR<=FDRcutoff)




VolcanoPlot<-ggplot(mainCorTable,aes(x=Pearsons,y=-log10(FDR+1e-16) ))+
geom_point()+
geom_hline(yintercept=-log10(FDRcutoff), color="red",size=1)+
geom_vline(xintercept=PearsonsCutoff, color="blue", size=1)+
geom_vline(xintercept=-PearsonsCutoff, color="blue", size=1)+
ylab("-log10( FDR )")+
xlab("Pearson's correlation")+
xlim(c(-1,1))

ggsave("VolcanoPlot.tiff", VolcanoPlot)



POEQTLVolcanoPlot<-ggplot(POEQTLmainCorTable,aes(x=Pearsons,y=-log10(FDR+1e-16) ))+
geom_point()+
geom_hline(yintercept=-log10(FDRcutoff), color="red",size=1)+
geom_vline(xintercept=PearsonsCutoff, color="blue", size=1)+
geom_vline(xintercept=-PearsonsCutoff, color="blue", size=1)+
ylab("-log10( FDR )")+
xlab("Pearson's correlation")+
xlim(c(-1,1))

ggsave("POEQTLVolcanoPlot.tiff", POEQTLVolcanoPlot)






CorHeatmap<-ggplot(mainCorTable,aes(x=Gene,y=Phenotype,fill=Pearsons))+
geom_tile()+
scale_fill_viridis(limits=c(-1,1))+
theme(axis.text.x=element_text(angle=90, hjust=1))

ggsave("CorHeatmap.tiff", CorHeatmap)



POEQTLCorHeatmap<-ggplot(POEQTLmainCorTable,aes(x=Gene,y=Phenotype,fill=Pearsons))+
geom_tile()+
scale_fill_viridis(limits=c(-1,1))+
theme(axis.text.x=element_text(angle=90, hjust=1))

ggsave("POEQTLCorHeatmap.tiff", POEQTLCorHeatmap)



unique(as.character(filtmainCorTable$Phenotype))


ListOfGenes<-unique(as.character(filtmainCorTable$Gene))


dim(EpiNet)

filtEpiNet<-subset(EpiNet,gA %in% ListOfGenes & gB %in% ListOfGenes )


dim(filtEpiNet)




write.table(filtmainCorTable, file= "filtmainCorTable.tsv", quote=FALSE, sep='\t', row.names=FALSE)







PhenNames<-c("Necweight","glu0","AUC.G.","TG_Serum","Cholesterol","GTTweight","averageFat","Insulin")


conversions<-list(
c("NecWeight","Go.dal","Ing","Kid_Fat","Mesen","Tot_Fat"),
c("T0_1","T0_2"),
c("AUC2","AUC1"),
c('TG'),
c("CHOL"),
c("GlucWt1","GlucWt2"),
c("NecWeight","Go.dal","Ing","Kid_Fat","Mesen","Tot_Fat"),
c("Insulin")
)



convfiltmainCorTable<-NULL

for(LINE in 1:nrow(filtmainCorTable)){

	PHENKEY<-as.character(filtmainCorTable[LINE,]$Phenotype)
	
	for(convID in conversions[PhenNames==PHENKEY][[1]]){
		#print(c(convID, PHENKEY,))
		
		convfiltmainCorTable<-rbind(cbind(filtmainCorTable[LINE,], convID), convfiltmainCorTable)


	}
		
}






HOLDER<-NULL


for(G in 1:length(ListOfGenes)){
	
	
	PHENDF<-subset(convfiltmainCorTable,Gene == ListOfGenes[G])
	
	for(COUNT in 1:nrow(PHENDF)){
		
		if( nrow(subset(filtEpiNet,gA == ListOfGenes[G])) >0 ){
		
			HOLDER<-rbind(cbind(PHENDF[COUNT,],subset(filtEpiNet,gA == ListOfGenes[G])), HOLDER)
		
		}
		
	}

}


for(G in 1:length(ListOfGenes)){
	
	
	PHENDF<-subset(convfiltmainCorTable,Gene == ListOfGenes[G])
	
	for(COUNT in 1:nrow(PHENDF)){
		
		if( nrow(subset(filtEpiNet,gB == ListOfGenes[G])) >0 ){
		
			HOLDER<-rbind(cbind(PHENDF[COUNT,],subset(filtEpiNet,gB == ListOfGenes[G])), HOLDER)
		
		}
		
	}

}



NetTable<-HOLDER[c(16,17,18,2,5)]

colnames(NetTable)<-c("Target Gene","ASE Gene","Pval","F1.Trait","Trait")

NetTable<-unique(NetTable)

write.table(NetTable, file= NETWORK_TABLE_FILE, quote=FALSE, sep='\t', row.names=FALSE)

