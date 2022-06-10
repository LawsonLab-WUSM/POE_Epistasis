

rm(list=ls())

#setwd("Box/LawsonLab/JuanMacias/scRNA_eWAT/")

library(mgcv)
library(dplyr)
library("Seurat")
library(patchwork)
library(cowplot)
library(ggplot2)
library(gridExtra)
library(viridis)
theme_set(theme_cowplot())
source("/Users/juanmacias/Box/LawsonLab/JuanMacias/Misc_Code/CS_Aesthetic.R")
library(tidyr)

data.seurat<-readRDS("data.seurat")


# SampleInformation<-read.delim("Sample_Information.tsv")
# 
# SampleInformation$Cohort<-paste(SampleInformation$Diet,SampleInformation$Sex,SampleInformation$Age,sep="")
# 
# SampleInformation$Index<-SampleInformation$IndexWithinBatch

### Based on nFeature < 3,000 threshold
new.cluster.ids <- c("MSC","Macrophage","VEC","VSMC","T cell","B cell","MSC.Diff")



ExpressionByCluster.Dataframe<-data.frame(
orig.ident=data.seurat[["orig.ident"]][,1],
Cluster=new.cluster.ids[data.seurat[["RNA_snn_res.0.15"]][,1]]
)




####### Finding cell type priors

ExpressionByCluster.Dataframe



ContextFreqTable<-table(ExpressionByCluster.Dataframe[,2])
ContextPriorTable<-ContextFreqTable/sum(ContextFreqTable)

ContextPriorTable




library(fitdistrplus)


new.cluster.ids

MSC<-subset(data.seurat, idents ="MSC")

Macrophage<-subset(data.seurat, idents = "Macrophage")

VEC<-subset(data.seurat, idents = "VEC")

VSMC<-subset(data.seurat, idents ="VSMC")

T.cell<-subset(data.seurat, idents = "T cell")

B.cell<-subset(data.seurat, idents ="B cell")

MSC.Diff<-subset(data.seurat, idents ="MSC.Diff")




MSC[["RNA"]]["Nnat"]



table(GetAssayData(MSC.Diff)[rownames(GetAssayData(MSC.Diff))=="Nnat",])
table(GetAssayData(MSC)[rownames(GetAssayData(MSC))=="Nnat",])



UMAPPlot(data.seurat)



FeaturePlot(data.seurat,features = "Nnat",order = TRUE)

FeaturePlot(data.seurat,features = "Adipoq",order = TRUE)



FeaturePlot(MSC.Diff,features = "Nnat",order = TRUE)

FeaturePlot(MSC.Diff,features = "Adipoq",order = TRUE)+FeaturePlot(MSC,features = "Adipoq",order = TRUE)




rownames(MSC[["RNA"]])

rownames(VEC[["RNA"]])





FindMaxandMin<-function(SEURAT.OBJECT,GENE,CELLTYPE.NAME){
  COUNTS<-as.numeric(SEURAT.OBJECT[["RNA"]][rownames(SEURAT.OBJECT[["RNA"]]) == GENE,])
  MIN<-min(COUNTS)
  MAX<-max(COUNTS)
  
  return(data.frame(MaxCount=MAX,MinCount=MIN,Gene=GENE,CellType=CELLTYPE.NAME))
}

FindMaxandMinForAllTypes<-function(SEURAT.OBJECT.LIST,GENE,CELLTYPE.NAME.LIST){
  
  lapply(1:length(SEURAT.OBJECT.LIST),function(TYPE) FindMaxandMin(SEURAT.OBJECT.LIST[[as.numeric(TYPE)]],GENE,CELLTYPE.NAME.LIST[TYPE]) ) %>% bind_rows()
  
}








EstimateDropoutAndFit<-function(GENE,SEURAT.OBJECT){
  COUNTS<-as.numeric(SEURAT.OBJECT[["RNA"]][rownames(SEURAT.OBJECT[["RNA"]]) == GENE,])
  COUNTS<-(COUNTS-min(COUNTS))/max(COUNTS)
  #Prior<-sum(COUNTS==0)/length(COUNTS)
  
  Conditionals<-seq(0,1,0.05)
  print(GENE)
  lapply(Conditionals, function(Pr.Dropout.Zero) runSingleFit(COUNTS, Pr.Dropout.Zero, GENE)) %>% bind_rows()
  
}



runSingleFit<-function(COUNTS, Pr.Dropout.Zero, GENE){
  
  ExpectedNumberOfDropouts<-ceiling(sum(COUNTS==0)*Pr.Dropout.Zero)
  NumberOfZeroCounts<-sum(COUNTS==0)
  
  
  COUNTS[COUNTS==0][sample(c(1:NumberOfZeroCounts),ExpectedNumberOfDropouts)]<-NA
  
  
  COUNTS<-COUNTS[!is.na(COUNTS)]
  
  FIT<-fitdist(COUNTS, distr = "beta", method = "mge")
  GOF<-gofstat(FIT)

  Shape1<-FIT$estimate[1]
  Shape2<-FIT$estimate[2]
  Ks<-GOF$ks
  
  OUTPUT<-data.frame(Shape1,Shape2,Ks,Pr.Dropout.Zero,GENE)
  rownames(OUTPUT)<-NULL
  return(OUTPUT)

}


findGenesMinNonZero<-function(COUNTS,GENE){
  PercentNonZero<-sum(as.numeric(COUNTS[["RNA"]][rownames(COUNTS[["RNA"]]) == GENE,])>0)/length(as.numeric(COUNTS[["RNA"]][rownames(COUNTS[["RNA"]]) == GENE,]))
  OUT<-data.frame(GENE,PercentNonZero)
  return(OUT)
}








# Gene.List<-rownames(Beta.1[["RNA"]])
# 
# NonZeroPercent.DF<-lapply(Gene.List, function(GENE) findGenesMinNonZero(Beta.1,GENE)) %>% bind_rows()
# 
# 
# Omit.List<-subset(NonZeroPercent.DF, PercentNonZero<0.05)$GENE
# Gene.List<-Gene.List[!(Gene.List %in% Omit.List)]
# 
# Fits.Out<-lapply(Gene.List, function(GENE) EstimateDropoutAndFit(GENE,Beta.1)) %>% bind_rows()

runFitForCellType<-function(SEURAT.OBJECT, MIN.NON.ZERO.FRACTION=0.05){
  
  print( c("Number of Cells: ",nrow(SEURAT.OBJECT[["RNA_snn_res.0.15"]])) )
  
  print("Extracting gene names...")
  GENE.LIST<-rownames(SEURAT.OBJECT[["RNA"]])
  print("Calculating Non zero fractions...")
  NonZeroPercent.DF<-lapply(GENE.LIST, function(GENE) findGenesMinNonZero(SEURAT.OBJECT,GENE)) %>% bind_rows()
  
  print("Deciding genes to omit...")
  OMIT.LIST<-subset(NonZeroPercent.DF, PercentNonZero<MIN.NON.ZERO.FRACTION)$GENE
  
  print(c("Total Genes: ",length(GENE.LIST)))
  print(c("Omitted Genes: ",length(OMIT.LIST)))
  GENE.LIST<-GENE.LIST[!(GENE.LIST %in% OMIT.LIST)]
  
  print("Running fitting...")
  OUT<-lapply(GENE.LIST, function(GENE) EstimateDropoutAndFit(GENE,SEURAT.OBJECT)) %>% bind_rows()
  print("Complete!")
  
  return(OUT)
}


#NonZeroPercent.DF


runFitForCellType<-function(SEURAT.OBJECT, NONZEROPERCENT.DF, MIN.NON.ZERO.FRACTION=0.05){
  
  print( c("Number of Cells: ",nrow(SEURAT.OBJECT[["RNA_snn_res.0.15"]])) )
  
  print("Extracting gene names...")
  GENE.LIST<-rownames(SEURAT.OBJECT[["RNA"]])
  print("Calculating Non zero fractions...")
  #NonZeroPercent.DF<-lapply(GENE.LIST, function(GENE) findGenesMinNonZero(SEURAT.OBJECT,GENE)) %>% bind_rows()
  
  print("Deciding genes to omit...")
  OMIT.LIST<-subset(NONZEROPERCENT.DF, PercentNonZero<MIN.NON.ZERO.FRACTION)$GENE
  
  print(c("Total Genes: ",length(GENE.LIST)))
  print(c("Omitted Genes: ",length(OMIT.LIST)))
  GENE.LIST<-GENE.LIST[!(GENE.LIST %in% OMIT.LIST)]
  
  print("Running fitting...")
  OUT<-lapply(GENE.LIST, function(GENE) EstimateDropoutAndFit(GENE,SEURAT.OBJECT)) %>% bind_rows()
  print("Complete!")
  
  return(OUT)
}


measureAllNonZeroFraction<-function(SEURAT.OBJECT){
  print("Extracting gene names...")
  GENE.LIST<-rownames(SEURAT.OBJECT[["RNA"]])
  print("Measuring fraction of non-zero cells...")
  NONZEROPERCENT.DF<-lapply(GENE.LIST, function(GENE) findGenesMinNonZero(SEURAT.OBJECT,GENE)) %>% bind_rows()
  print("Complete!")
  return(NONZEROPERCENT.DF)
}





table(as.numeric(MSC[["RNA"]]["Nnat"]))

table(as.numeric(MSC.Diff[["RNA"]]["Nnat"]))



#subset(NonZeroFraction.MSC,GENE=="Adipoq")


# NonZeroFraction.MSC<-measureAllNonZeroFraction(MSC)
# Fits.MSC<-runFitForCellType(MSC,NonZeroFraction.MSC)
# write.table(Fits.MSC,file = "Fits.MSC.txt",quote = FALSE,col.names = TRUE, row.names = FALSE, sep ="\t")
# write.table(NonZeroFraction.MSC,file = "NonZeroFraction.MSC.txt",quote = FALSE,col.names = TRUE, row.names = FALSE, sep ="\t")
# 
# NonZeroFraction.Macrophage<-measureAllNonZeroFraction(Macrophage)
# Fits.Macrophage<-runFitForCellType(Macrophage,NonZeroFraction.Macrophage)
# write.table(Fits.Macrophage,file = "Fits.Macrophage.txt",quote = FALSE,col.names = TRUE, row.names = FALSE, sep ="\t")
# write.table(NonZeroFraction.Macrophage,file = "NonZeroFraction.Macrophage.txt",quote = FALSE,col.names = TRUE, row.names = FALSE, sep ="\t")
# 
# NonZeroFraction.VEC<-measureAllNonZeroFraction(VEC)
# Fits.VEC<-runFitForCellType(VEC,NonZeroFraction.VEC)
# write.table(Fits.VEC,file = "Fits.VEC.txt",quote = FALSE,col.names = TRUE, row.names = FALSE, sep ="\t")
# write.table(NonZeroFraction.VEC,file = "NonZeroFraction.VEC.txt",quote = FALSE,col.names = TRUE, row.names = FALSE, sep ="\t")
# 
# NonZeroFraction.VSMC<-measureAllNonZeroFraction(VSMC)
# Fits.VSMC<-runFitForCellType(VSMC,NonZeroFraction.VSMC)
# write.table(Fits.VSMC,file = "Fits.VSMC.txt",quote = FALSE,col.names = TRUE, row.names = FALSE, sep ="\t")
# write.table(NonZeroFraction.VSMC,file = "NonZeroFraction.VSMC.txt",quote = FALSE,col.names = TRUE, row.names = FALSE, sep ="\t")
# 
# NonZeroFraction.T.cell<-measureAllNonZeroFraction(T.cell)
# Fits.T.cell<-runFitForCellType(T.cell,NonZeroFraction.T.cell)
# write.table(Fits.T.cell,file = "Fits.T.cell.txt",quote = FALSE,col.names = TRUE, row.names = FALSE, sep ="\t")
# write.table(NonZeroFraction.T.cell,file = "NonZeroFraction.T.cell.txt",quote = FALSE,col.names = TRUE, row.names = FALSE, sep ="\t")
# 
# NonZeroFraction.B.cell<-measureAllNonZeroFraction(B.cell)
# Fits.B.cell<-runFitForCellType(B.cell,NonZeroFraction.B.cell)
# write.table(Fits.B.cell,file = "Fits.B.cell.txt",quote = FALSE,col.names = TRUE, row.names = FALSE, sep ="\t")
# write.table(NonZeroFraction.B.cell,file = "NonZeroFraction.B.cell.txt",quote = FALSE,col.names = TRUE, row.names = FALSE, sep ="\t")
# 
# NonZeroFraction.MSC.Diff<-measureAllNonZeroFraction(MSC.Diff)
# Fits.MSC.Diff<-runFitForCellType(MSC.Diff,NonZeroFraction.MSC.Diff)
# write.table(Fits.MSC.Diff,file = "Fits.MSC.Diff.txt",quote = FALSE,col.names = TRUE, row.names = FALSE, sep ="\t")
# write.table(NonZeroFraction.MSC.Diff,file = "NonZeroFraction.MSC.Diff.txt",quote = FALSE,col.names = TRUE, row.names = FALSE, sep ="\t")


new.cluster.ids

########


NonZeroFraction.MSC<-read.delim("NonZeroFraction.MSC.txt")
Fits.MSC<-read.delim("Fits.MSC.txt")

NonZeroFraction.Macrophage<-read.delim("NonZeroFraction.Macrophage.txt")
Fits.Macrophage<-read.delim("Fits.Macrophage.txt")

NonZeroFraction.VEC<-read.delim("NonZeroFraction.VEC.txt")
Fits.VEC<-read.delim("Fits.VEC.txt")

NonZeroFraction.VSMC<-read.delim("NonZeroFraction.VSMC.txt")
Fits.VSMC<-read.delim("Fits.VSMC.txt")

NonZeroFraction.T.cell<-read.delim("NonZeroFraction.T.cell.txt")
Fits.T.cell<-read.delim("Fits.T.cell.txt")

NonZeroFraction.B.cell<-read.delim("NonZeroFraction.B.cell.txt")
Fits.B.cell<-read.delim("Fits.B.cell.txt")

NonZeroFraction.MSC.Diff<-read.delim("NonZeroFraction.MSC.Diff.txt")
Fits.MSC.Diff<-read.delim("Fits.MSC.Diff.txt")


dim(NonZeroFraction.MSC)

NonZeroFraction.All<-rbind(
  NonZeroFraction.MSC,
  NonZeroFraction.Macrophage,
  NonZeroFraction.VEC,
  NonZeroFraction.VSMC,
  NonZeroFraction.T.cell,
  NonZeroFraction.B.cell,
  NonZeroFraction.MSC.Diff
  )

dim(NonZeroFraction.All)


sum(table(NonZeroFraction.All$GENE) !=7)

NonZeroFraction.All$CellType<-c(
  rep("MSC",nrow(NonZeroFraction.MSC)),
  rep("Macrophage",nrow(NonZeroFraction.Macrophage)),
  rep("VEC",nrow(NonZeroFraction.VEC)),
  rep("VSMC",nrow(NonZeroFraction.VSMC)),
  rep("T Cell",nrow(NonZeroFraction.T.cell)),
  rep("B Cell",nrow(NonZeroFraction.B.cell)),
  rep("Diff MSC",nrow(NonZeroFraction.MSC.Diff))
)

dim(NonZeroFraction.All)


NonZeroFraction.All.Wide<-spread(NonZeroFraction.All,CellType,PercentNonZero)
NonZeroFraction.All.Matrix<-as.matrix(NonZeroFraction.All.Wide[,-1])
rownames(NonZeroFraction.All.Matrix)<-NonZeroFraction.All.Wide$GENE


NonZeroFraction.Gene.dendro <- as.dendrogram(hclust(d = dist(x = NonZeroFraction.All.Matrix)))
NonZeroFraction.Gene.order <- order.dendrogram(NonZeroFraction.Gene.dendro)
NonZeroFraction.All$GENE<-factor(NonZeroFraction.All$GENE,levels = rownames(NonZeroFraction.All.Matrix)[NonZeroFraction.Gene.order])

NonZeroFraction.Type.dendro <- as.dendrogram(hclust(d = dist(x = t(NonZeroFraction.All.Matrix) )))
NonZeroFraction.Type.order <- order.dendrogram(NonZeroFraction.Type.dendro)
NonZeroFraction.All$CellType<-factor(NonZeroFraction.All$CellType,levels = colnames(NonZeroFraction.All.Matrix)[NonZeroFraction.Type.order])

NonZeroFraction.Heatmap<-ggplot(NonZeroFraction.All,aes(x=GENE,y=CellType,fill=PercentNonZero))+
  geom_tile()+
  theme(axis.text.x = element_blank())+
  scale_fill_viridis()+
  CS.THEME

ggsave(filename = "NonZeroFraction.Heatmap.png",plot = NonZeroFraction.Heatmap,units = "in",width = 24, height = 6)


######## NOTE: In cases where several parameter combinations fit equally well, we choose the set with the smallest Pr.Dropout.Zero, which should be the first row

# Fits.Optimal.MSC<-lapply(unique(Fits.MSC$GENE), function(Gene) subset(subset(Fits.MSC, GENE==Gene),Ks==min(subset(Fits.MSC, GENE==Gene)$Ks, na.rm = TRUE))[1,] ) %>% bind_rows()
# write.table(Fits.Optimal.MSC,file = "Fits.Optimal.MSC.txt",quote = FALSE,col.names = TRUE, row.names = FALSE, sep ="\t")
# 
# Fits.Optimal.Macrophage<-lapply(unique(Fits.Macrophage$GENE), function(Gene) subset(subset(Fits.Macrophage, GENE==Gene),Ks==min(subset(Fits.Macrophage, GENE==Gene)$Ks, na.rm = TRUE))[1,] ) %>% bind_rows()
# write.table(Fits.Optimal.Macrophage,file = "Fits.Optimal.Macrophage.txt",quote = FALSE,col.names = TRUE, row.names = FALSE, sep ="\t")
# 
# Fits.Optimal.VEC<-lapply(unique(Fits.VEC$GENE), function(Gene) subset(subset(Fits.VEC, GENE==Gene),Ks==min(subset(Fits.VEC, GENE==Gene)$Ks, na.rm = TRUE))[1,] ) %>% bind_rows()
# write.table(Fits.Optimal.VEC,file = "Fits.Optimal.VEC.txt",quote = FALSE,col.names = TRUE, row.names = FALSE, sep ="\t")
# 
# Fits.Optimal.VSMC<-lapply(unique(Fits.VSMC$GENE), function(Gene) subset(subset(Fits.VSMC, GENE==Gene),Ks==min(subset(Fits.VSMC, GENE==Gene)$Ks, na.rm = TRUE))[1,] ) %>% bind_rows()
# write.table(Fits.Optimal.VSMC,file = "Fits.Optimal.VSMC.txt",quote = FALSE,col.names = TRUE, row.names = FALSE, sep ="\t")
# 
# Fits.Optimal.T.cell<-lapply(unique(Fits.T.cell$GENE), function(Gene) subset(subset(Fits.T.cell, GENE==Gene),Ks==min(subset(Fits.T.cell, GENE==Gene)$Ks, na.rm = TRUE))[1,] ) %>% bind_rows()
# write.table(Fits.Optimal.T.cell,file = "Fits.Optimal.T.cell.txt",quote = FALSE,col.names = TRUE, row.names = FALSE, sep ="\t")
# 
# Fits.Optimal.B.cell<-lapply(unique(Fits.B.cell$GENE), function(Gene) subset(subset(Fits.B.cell, GENE==Gene),Ks==min(subset(Fits.B.cell, GENE==Gene)$Ks, na.rm = TRUE))[1,] ) %>% bind_rows()
# write.table(Fits.Optimal.B.cell,file = "Fits.Optimal.B.cell.txt",quote = FALSE,col.names = TRUE, row.names = FALSE, sep ="\t")
# 
# Fits.Optimal.MSC.Diff<-lapply(unique(Fits.MSC.Diff$GENE), function(Gene) subset(subset(Fits.MSC.Diff, GENE==Gene),Ks==min(subset(Fits.MSC.Diff, GENE==Gene)$Ks, na.rm = TRUE))[1,] ) %>% bind_rows()
# write.table(Fits.Optimal.MSC.Diff,file = "Fits.Optimal.MSC.Diff.txt",quote = FALSE,col.names = TRUE, row.names = FALSE, sep ="\t")





Fits.Optimal.MSC<-read.delim("Fits.Optimal.MSC.txt")
Fits.Optimal.Macrophage<-read.delim("Fits.Optimal.Macrophage.txt")
Fits.Optimal.VEC<-read.delim("Fits.Optimal.VEC.txt")
Fits.Optimal.VSMC<-read.delim("Fits.Optimal.VSMC.txt")
Fits.Optimal.T.cell<-read.delim("Fits.Optimal.T.cell.txt")
Fits.Optimal.B.cell<-read.delim("Fits.Optimal.B.cell.txt")
Fits.Optimal.MSC.Diff<-read.delim("Fits.Optimal.MSC.Diff.txt")


### Here we are integrating Pr(0) data with optimal fit data


NonZeroFraction.MSC$Pr.Zero<-1-NonZeroFraction.MSC$PercentNonZero
NonZeroFraction.Macrophage$Pr.Zero<-1-NonZeroFraction.Macrophage$PercentNonZero
NonZeroFraction.VEC$Pr.Zero<-1-NonZeroFraction.VEC$PercentNonZero
NonZeroFraction.VSMC$Pr.Zero<-1-NonZeroFraction.VSMC$PercentNonZero
NonZeroFraction.T.cell$Pr.Zero<-1-NonZeroFraction.T.cell$PercentNonZero
NonZeroFraction.B.cell$Pr.Zero<-1-NonZeroFraction.B.cell$PercentNonZero
NonZeroFraction.MSC.Diff$Pr.Zero<-1-NonZeroFraction.MSC.Diff$PercentNonZero





Fits.Optimal.MSC$Pr.Zero<-NA
Fits.Optimal.Macrophage$Pr.Zero<-NA
Fits.Optimal.VEC$Pr.Zero<-NA
Fits.Optimal.VSMC$Pr.Zero<-NA
Fits.Optimal.T.cell$Pr.Zero<-NA
Fits.Optimal.B.cell$Pr.Zero<-NA
Fits.Optimal.MSC.Diff$Pr.Zero<-NA




Fits.Optimal.MSC[match(subset(NonZeroFraction.MSC,GENE %in% Fits.Optimal.MSC$GENE)$GENE, Fits.Optimal.MSC$GENE),]$Pr.Zero<-subset(NonZeroFraction.MSC,GENE %in% Fits.Optimal.MSC$GENE)$Pr.Zero
Fits.Optimal.Macrophage[match(subset(NonZeroFraction.Macrophage,GENE %in% Fits.Optimal.Macrophage$GENE)$GENE, Fits.Optimal.Macrophage$GENE),]$Pr.Zero<-subset(NonZeroFraction.Macrophage,GENE %in% Fits.Optimal.Macrophage$GENE)$Pr.Zero
Fits.Optimal.VEC[match(subset(NonZeroFraction.VEC,GENE %in% Fits.Optimal.VEC$GENE)$GENE, Fits.Optimal.VEC$GENE),]$Pr.Zero<-subset(NonZeroFraction.VEC,GENE %in% Fits.Optimal.VEC$GENE)$Pr.Zero
Fits.Optimal.VSMC[match(subset(NonZeroFraction.VSMC,GENE %in% Fits.Optimal.VSMC$GENE)$GENE, Fits.Optimal.VSMC$GENE),]$Pr.Zero<-subset(NonZeroFraction.VSMC,GENE %in% Fits.Optimal.VSMC$GENE)$Pr.Zero
Fits.Optimal.T.cell[match(subset(NonZeroFraction.T.cell,GENE %in% Fits.Optimal.T.cell$GENE)$GENE, Fits.Optimal.T.cell$GENE),]$Pr.Zero<-subset(NonZeroFraction.T.cell,GENE %in% Fits.Optimal.T.cell$GENE)$Pr.Zero
Fits.Optimal.B.cell[match(subset(NonZeroFraction.B.cell,GENE %in% Fits.Optimal.B.cell$GENE)$GENE, Fits.Optimal.B.cell$GENE),]$Pr.Zero<-subset(NonZeroFraction.B.cell,GENE %in% Fits.Optimal.B.cell$GENE)$Pr.Zero
Fits.Optimal.MSC.Diff[match(subset(NonZeroFraction.MSC.Diff,GENE %in% Fits.Optimal.MSC.Diff$GENE)$GENE, Fits.Optimal.MSC.Diff$GENE),]$Pr.Zero<-subset(NonZeroFraction.MSC.Diff,GENE %in% Fits.Optimal.MSC.Diff$GENE)$Pr.Zero





Fits.Optimal.MSC$Pr.Dropout<-Fits.Optimal.MSC$Pr.Zero*Fits.Optimal.MSC$Pr.Dropout.Zero
Fits.Optimal.Macrophage$Pr.Dropout<-Fits.Optimal.Macrophage$Pr.Zero*Fits.Optimal.Macrophage$Pr.Dropout.Zero
Fits.Optimal.VEC$Pr.Dropout<-Fits.Optimal.VEC$Pr.Zero*Fits.Optimal.VEC$Pr.Dropout.Zero
Fits.Optimal.VSMC$Pr.Dropout<-Fits.Optimal.VSMC$Pr.Zero*Fits.Optimal.VSMC$Pr.Dropout.Zero
Fits.Optimal.T.cell$Pr.Dropout<-Fits.Optimal.T.cell$Pr.Zero*Fits.Optimal.T.cell$Pr.Dropout.Zero
Fits.Optimal.B.cell$Pr.Dropout<-Fits.Optimal.B.cell$Pr.Zero*Fits.Optimal.B.cell$Pr.Dropout.Zero
Fits.Optimal.MSC.Diff$Pr.Dropout<-Fits.Optimal.MSC.Diff$Pr.Zero*Fits.Optimal.MSC.Diff$Pr.Dropout.Zero






Fits.Optimal.MSC$CellType<-"MSC"
Fits.Optimal.Macrophage$CellType<-"Macrophage"
Fits.Optimal.VEC$CellType<-"VEC"
Fits.Optimal.VSMC$CellType<-"VSMC"
Fits.Optimal.T.cell$CellType<-"T cell"
Fits.Optimal.B.cell$CellType<-"B cell"
Fits.Optimal.MSC.Diff$CellType<-"MSC.Diff"


Fits.Optimal.All<-rbind(
  Fits.Optimal.MSC,
  Fits.Optimal.Macrophage,
  Fits.Optimal.VEC,
  Fits.Optimal.VSMC,
  Fits.Optimal.T.cell,
  Fits.Optimal.B.cell,
  Fits.Optimal.MSC.Diff
  )

Fits.Optimal.All$CellType<-factor(Fits.Optimal.All$CellType,names(ContextPriorTable))

Fits.Optimal.All$ExpectedValue<-Fits.Optimal.All$Shape1/(Fits.Optimal.All$Shape1+Fits.Optimal.All$Shape2)


MeanError.Plot.MSC<-ggplot(subset(Fits.Optimal.All,CellType=="MSC"),aes(x=Ks,y=ExpectedValue,color=Pr.Zero))+
  geom_point(size=0.1)+
  scale_y_continuous(limits = c(0,1))+
  scale_x_continuous(limits = c(0,0.15))+
  ggtitle("Mesenchymal Stem Cell")+
  scale_color_viridis(end = 0.9,option = "plasma")+
  ylab("Expected value")+
  xlab("K-S Statistic")+
  CS.THEME

MeanError.Legend<-get_legend(MeanError.Plot.MSC)

MeanError.Plot.MSC<-MeanError.Plot.MSC+
  theme(legend.position = "none")

MeanError.Plot.Macrophage<-ggplot(subset(Fits.Optimal.All,CellType=="Macrophage"),aes(x=Ks,y=ExpectedValue,color=Pr.Zero))+
  geom_point(size=0.1)+
  scale_y_continuous(limits = c(0,1))+
  scale_x_continuous(limits = c(0,0.15))+
  ggtitle("Macrophage")+
  scale_color_viridis(end = 0.9,option = "plasma")+
  ylab("Expected value")+
  xlab("K-S Statistic")+
  CS.THEME+
  theme(legend.position = "none")

MeanError.Plot.VEC<-ggplot(subset(Fits.Optimal.All,CellType=="VEC"),aes(x=Ks,y=ExpectedValue,color=Pr.Zero))+
  geom_point(size=0.1)+
  scale_y_continuous(limits = c(0,1))+
  scale_x_continuous(limits = c(0,0.15))+
  ggtitle("Vascular Endothelial Cell")+
  scale_color_viridis(end = 0.9,option = "plasma")+
  ylab("Expected value")+
  xlab("K-S Statistic")+
  CS.THEME+
  theme(legend.position = "none")

MeanError.Plot.VSMC<-ggplot(subset(Fits.Optimal.All,CellType=="VSMC"),aes(x=Ks,y=ExpectedValue,color=Pr.Zero))+
  geom_point(size=0.1)+
  scale_y_continuous(limits = c(0,1))+
  scale_x_continuous(limits = c(0,0.15))+
  ggtitle("Vascular Smooth Muscle Cell")+
  scale_color_viridis(end = 0.9,option = "plasma")+
  ylab("Expected value")+
  xlab("K-S Statistic")+
  CS.THEME+
  theme(legend.position = "none")

MeanError.Plot.T.Cell<-ggplot(subset(Fits.Optimal.All,CellType=="T cell"),aes(x=Ks,y=ExpectedValue,color=Pr.Zero))+
  geom_point(size=0.1)+
  scale_y_continuous(limits = c(0,1))+
  scale_x_continuous(limits = c(0,0.15))+
  ggtitle("T cell")+
  scale_color_viridis(end = 0.9,option = "plasma")+
  ylab("Expected value")+
  xlab("K-S Statistic")+
  CS.THEME+
  theme(legend.position = "none")

MeanError.Plot.B.Cell<-ggplot(subset(Fits.Optimal.All,CellType=="B cell"),aes(x=Ks,y=ExpectedValue,color=Pr.Zero))+
  geom_point(size=0.1)+
  scale_y_continuous(limits = c(0,1))+
  scale_x_continuous(limits = c(0,0.15))+
  ggtitle("B cell")+
  scale_color_viridis(end = 0.9,option = "plasma")+
  ylab("Expected value")+
  xlab("K-S Statistic")+
  CS.THEME+
  theme(legend.position = "none")

MeanError.Plot.MSC.Diff<-ggplot(subset(Fits.Optimal.All,CellType=="MSC.Diff"),aes(x=Ks,y=ExpectedValue,color=Pr.Zero))+
  geom_point(size=0.1)+
  scale_y_continuous(limits = c(0,1))+
  scale_x_continuous(limits = c(0,0.15))+
  ggtitle("Differentiating Mesenchymal Stem Cell")+
  scale_color_viridis(end = 0.9,option = "plasma")+
  ylab("Expected value")+
  xlab("K-S Statistic")+
  CS.THEME+
  theme(legend.position = "none")


MeanError.Panel<-plot_grid(
  MeanError.Plot.MSC,
  MeanError.Plot.MSC.Diff,
  MeanError.Plot.Macrophage,
  MeanError.Plot.VEC,
  MeanError.Plot.VSMC,
  MeanError.Plot.T.Cell,
  MeanError.Plot.B.Cell,
  MeanError.Legend,
  labels = c("A","B","C","D","E",'F',"G")
)

ggsave2(filename = "MeanError.Panel.png",plot = MeanError.Panel,units = "in",width = 16,height = 9)




Beta.Hurdle.Statistic_Pr.Dropout.Zero<-ggplot(Fits.Optimal.All,aes(x=Pr.Dropout.Zero))+
  geom_histogram()+
  xlab("Pr( Dropout | Zero )")

Beta.Hurdle.Statistic_Ks<-ggplot(Fits.Optimal.All,aes(x=Ks))+
  geom_histogram()+
  xlab("K-S Statistic")

Beta.Hurdle.Statistic_Pr.Zero<-ggplot(Fits.Optimal.All,aes(x=Pr.Zero))+
  geom_histogram()+
  xlab("Pr( Zero )")

Beta.Hurdle.Statistic_ExpectedValue<-ggplot(Fits.Optimal.All,aes(x=ExpectedValue))+
  geom_histogram()+
  xlab("Expected Value")



Beta.Hurdle.Statistic.Panel<-plot_grid(Beta.Hurdle.Statistic_Pr.Zero,Beta.Hurdle.Statistic_Pr.Dropout.Zero,Beta.Hurdle.Statistic_Ks,Beta.Hurdle.Statistic_ExpectedValue,labels = "AUTO")


ggsave2(filename = "Beta.Hurdle.Statistic.Panel.png",plot = Beta.Hurdle.Statistic.Panel,units = "in",width = 7,height = 6)


### MCSM

Seurat.Object.List.All<-c(MSC,Macrophage,VEC,VSMC,T.cell,B.cell,MSC.Diff)



calculateExpectedContributions<-function(FITS.OPTIMAL,SEURAT.OBJECT.LIST,GENE.QUERY,CELLTYPE.PR.NAMES,CELLTYPE.PR.VECTOR){
  TEST<-subset(FITS.OPTIMAL,GENE==GENE.QUERY)
  
  TEST.Scales<-FindMaxandMinForAllTypes(SEURAT.OBJECT.LIST,GENE.QUERY,CELLTYPE.PR.NAMES)
  
  
  CellTypeProbabilities<-CELLTYPE.PR.VECTOR[match(TEST$CellType,CELLTYPE.PR.NAMES)]
  CellTypeNonDropoutProbabilities<-(1-TEST$Pr.Dropout)
  CellTypeExpectedValues<-TEST$ExpectedValue
  CellTypeMax<-TEST.Scales$MaxCount[match(TEST$CellType,TEST.Scales$CellType)]
  CellTypeMin<-TEST.Scales$MinCount[match(TEST$CellType,TEST.Scales$CellType)]
  
  CONTRIBUTIONS<-CellTypeProbabilities*CellTypeNonDropoutProbabilities*(CellTypeExpectedValues*CellTypeMax-CellTypeMin)
  
  CONTRIBUTIONS.DF<-data.frame(
    CellType=CELLTYPE.PR.NAMES,
    Contribution=NA
  )
  
  CONTRIBUTIONS.DF[match(names(CONTRIBUTIONS),CONTRIBUTIONS.DF$CellType),]$Contribution<-CONTRIBUTIONS
  
  CONTRIBUTIONS.DF$Contribution[is.na(CONTRIBUTIONS.DF$Contribution)]<-0
  
  CONTRIBUTIONS.DF$FractionOfTotal<-CONTRIBUTIONS.DF$Contribution/sum(CONTRIBUTIONS.DF$Contribution)
  
  CONTRIBUTIONS.DF$Gene<-GENE.QUERY
  
  CONTRIBUTIONS.DF
}



Contributions<-lapply(
  unique(as.character(Fits.Optimal.All$GENE)),
  function(GENEOFINTEREST) calculateExpectedContributions(
    Fits.Optimal.All,
    Seurat.Object.List.All,
    GENEOFINTEREST,
    names(ContextPriorTable),
    ContextPriorTable
  )
) %>% bind_rows()



write.table(Contributions,"Contributions.txt",row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")



CellTypePriorsFromLit<-ContextPriorTable

CellTypePriorsFromLit[1]<-0.6
CellTypePriorsFromLit[2]<-4.4
CellTypePriorsFromLit[3]<-15.7
CellTypePriorsFromLit[4]<-66.4
CellTypePriorsFromLit[5]<-2.2
CellTypePriorsFromLit[6]<-1.8
CellTypePriorsFromLit[7]<-8.9

CellTypePriorsFromLit<-CellTypePriorsFromLit/100



Contributions.Lit<-lapply(
  unique(as.character(Fits.Optimal.All$GENE)),
  function(GENEOFINTEREST) calculateExpectedContributions(
    Fits.Optimal.All,
    Seurat.Object.List.All,
    GENEOFINTEREST,
    names(ContextPriorTable),
    CellTypePriorsFromLit
  )
) %>% bind_rows()

write.table(Contributions.Lit,"Contributions.Literature.txt",row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")




library(tidyr)
Contributions.Lit<-read.delim("Contributions.Literature.txt")
Contributions.Lit.Wide<-spread(Contributions.Lit[,c(1,3,4)],CellType,FractionOfTotal)

write.table(Contributions.Lit.Wide,"Contributions.Lit.Wide.txt",row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")


Contributions<-read.delim("Contributions.txt")
Contributions.Wide<-spread(Contributions[,c(1,3,4)],CellType,FractionOfTotal)


write.table(Contributions.Wide,"Contributions.Wide.txt",row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")


FractionOfTotal.DistPlot<-ggplot(Contributions,aes(x=FractionOfTotal))+
  geom_histogram()+
  facet_grid(CellType~.)+
  CS.THEME+
  xlab("Relative contribution to expression")

ggsave2(filename = "FractionOfTotal.DistPlot.png",plot = FractionOfTotal.DistPlot,units = "in",width = 4,height = 9)


GOI<-c("Sfrp4","Reep6","Pde2a","Myl7","Ly6d","Csf3r","Cd28")
ExamplesUMAP<-FeaturePlot(data.seurat, features = GOI, order=TRUE, pt.size = 0.5,ncol = 4)+UMAPPlot(data.seurat)+CS.THEME


ggsave2(filename = "ExamplesUMAP.png",plot = ExamplesUMAP,units = "in",width = 12,height = 5)





