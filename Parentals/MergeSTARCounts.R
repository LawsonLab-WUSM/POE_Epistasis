
rm(list=ls())
set.seed(0)

require("edgeR")
library(ggplot2)
library(cowplot)
source("/Users/juanmacias/Box/LawsonLab/JuanMacias/Misc_Code/CS_Aesthetic.R")

SM<-read.table("SM_BAT_WAT_GeneCounts.tsv",header=TRUE)
SM.summary<-SM[1:4,2:ncol(SM)]
SM<-SM[5:nrow(SM),]

LG<-read.table("LG_BAT_WAT_GeneCounts.tsv",header=TRUE)
LG.summary<-LG[1:4,2:ncol(LG)]
LG<-LG[5:nrow(LG),]


Recip<-read.table("WAT_MergedCounts.tsv",header=TRUE)
Recip.summary<-Recip[1:4,2:ncol(Recip)]
Recip<-Recip[5:nrow(Recip),]

SM.Mixed<-read.table("SM_mixed_MergedCounts.tsv",header=TRUE)
SM.Mixed.summary<-SM.Mixed[1:4,2:ncol(SM.Mixed)]
SM.Mixed<-SM.Mixed[5:nrow(SM.Mixed),]

LG.Mixed<-read.table("LG_mixed_MergedCounts.tsv",header=TRUE)
LG.Mixed.summary<-LG.Mixed[1:4,2:ncol(LG.Mixed)]
LG.Mixed<-LG.Mixed[5:nrow(LG.Mixed),]





All.Summaries<-cbind(SM.summary,LG.summary,Recip.summary,SM.Mixed.summary, LG.Mixed.summary)
Run<-c(rep("SM",ncol(SM.summary)),rep("LG",ncol(LG.summary)),rep("Recip",ncol(Recip.summary)),rep("Mixed.SM",ncol(SM.Mixed.summary)),rep("Mixed.LG",ncol(LG.Mixed.summary)))

### N_unmapped, N_multimapping, N_noFeature, N_ambiguous

plot(x=t(All.Summaries[1,]),y=t(All.Summaries[2,]),pch=16,cex=0.5,col=c("red","darkgreen","green","orange","blue")[as.numeric(factor(Run))])

plot(x=t(All.Summaries[1,]),y=t(All.Summaries[4,]),pch=16,cex=0.5,col=c("red","darkgreen","green","orange","blue")[as.numeric(factor(Run))])


SM.Genes<-as.character(SM$Gene)
LG.Genes<-as.character(LG$Gene)
Recip.Genes<-as.character(Recip$Gene)
SM.Mixed.Genes<-as.character(SM.Mixed$Gene)
LG.Mixed.Genes<-as.character(LG.Mixed$Gene)



sum(!(SM.Mixed.Genes %in% LG.Mixed.Genes))
sum(!(Recip.Genes %in% LG.Mixed.Genes))
sum(!(Recip.Genes %in% SM.Mixed.Genes))

sum(!(LG.Mixed.Genes %in% SM.Mixed.Genes))
sum(!(LG.Mixed.Genes %in% Recip.Genes))
sum(!(SM.Mixed.Genes %in% Recip.Genes))

length(LG.Mixed.Genes)
length(SM.Mixed.Genes)
length(Recip.Genes)

sum(!(SM.Mixed.Genes == LG.Mixed.Genes))
sum(!(SM.Mixed.Genes ==Recip.Genes))


####### Based on the above I should be able to collapse down the above into one DF easily


DF<-cbind(SM.Mixed[,2:ncol(SM.Mixed)], LG.Mixed[,2:ncol(LG.Mixed)],Recip[,2:ncol(Recip)])

DF.Mat<-as.matrix(DF,ncol=ncol(DF))
rownames(DF.Mat)<-SM.Mixed.Genes





Sample.Info <-read.delim("All_BAT_WAT_SampleInfo.txt")




#Order Sample Info
Sample.Info<-Sample.Info[match(colnames(DF.Mat),as.character(Sample.Info$GTAC.Tag)),]

## WAT only
DF.Mat<-DF.Mat[,Sample.Info$Tissue=="r"]
Sample.Info<-Sample.Info[Sample.Info$Tissue=="r",]


fERdf<-DGEList(counts=DF.Mat)

#Filtering
minCPM<-10
minSample<-40

keep<-rowSums(fERdf$counts>=minCPM) >=minSample
fERdf<-fERdf[keep,keep.lib.sizes=FALSE]

#Normalization
fERdf<-calcNormFactors(fERdf,method="TMM")




LibrarySummaries<-data.frame(
  Sample=rownames(fERdf$samples),
  LibrarySize=fERdf$samples$lib.size,
  NormFactor=fERdf$samples$norm.factors,
  Genotype=Sample.Info$Genotype[match(rownames(fERdf$samples),Sample.Info$GTAC.Tag)]
)







LibraryMetricsPlot<-ggplot(LibrarySummaries,aes(x=LibrarySize/1000000,y=NormFactor,color=Genotype))+
  geom_point(size=2)+
  scale_color_viridis(discrete = TRUE,end=0.95,option = "plasma")+
  CS.THEME+
  xlab("Library size (M)")+
  ylab("Normalization factor")+
  theme(legend.position = "bottom")



ggsave(filename = "TotalExpression_LibraryMetricsPlot.png",plot = LibraryMetricsPlot,width = 3.25,height = 3.5)


SEX<-Sample.Info$Sex
DIET<-Sample.Info$Diet
AGE<-Sample.Info$Age
GENOTYPE<-Sample.Info$Genotype

design <- model.matrix(~ SEX*DIET*GENOTYPE+AGE)

fERdf<-estimateDisp(fERdf,design)

fit<-glmFit(fERdf, design)

#Save Counts
Counts<-cpm(fERdf, normalized.lib.sizes=TRUE)



write.table(Counts,"Normalized_Counts.tsv",sep="\t",row.names=TRUE,col.names=TRUE,quote=FALSE)



