
rm(list=ls())
library(mgcv)
library(dplyr)
library("Seurat")
library(patchwork)
library(cowplot)
library(ggplot2)
library(gridExtra)
library(viridis)
library(DoubletDecon)
theme_set(theme_cowplot())



theme_set(theme(
		panel.background = element_rect(fill="white"),
		panel.grid.major=element_line(color="gray90", linetype="solid"), #adds major grid lines
		panel.grid.minor=element_line(color="gray95", linetype="solid", size=0.65), #adds minor grid lines	
		panel.border=element_rect(fill=NA, color="black", size=1.5, linetype="solid"), #draws a black border around plot
		axis.line=element_blank()
	))

theme_set(theme(text=element_text(colour = "black", family= "Arial")))




filterMatrixDirectory <-"eWAT_All_filtered_feature_bc_matrix"

data <- Read10X(data.dir = filterMatrixDirectory)

data.seurat<-CreateSeuratObject(counts = data, min.cells = 200, min.features = 500)
data.seurat






minFeatures<-500
maxFeatures<-3000
minCount<-1000
maxCount<-30000


QC.Dataframe<-data.frame(
orig.ident=data.seurat[["orig.ident"]],
nFeature_RNA=data.seurat[["nFeature_RNA"]],
nCount_RNA=data.seurat[["nCount_RNA"]]
)




nFeatureDistribution.Plot<-ggplot(QC.Dataframe,aes(x= nFeature_RNA))+
geom_histogram(bins=70)+
geom_vline(xintercept= minFeatures,color="red",size=1,alpha=0.5)+
geom_vline(xintercept= maxFeatures,color="red",size=1,alpha=0.5)

ggsave2("nFeatureDistribution.Plot.png", nFeatureDistribution.Plot,units="in",width=6,height=3)




nCount.Plot<-ggplot(QC.Dataframe,aes(x= nCount_RNA))+
geom_histogram(bins=70)+
geom_vline(xintercept= minCount,color="red",size=1,alpha=0.5)+
geom_vline(xintercept= maxCount,color="red",size=1,alpha=0.5)+
scale_x_continuous(trans='log2')

ggsave2("nCount.Plot.png", nCount.Plot,units="in",width=6,height=3)




Principal_QC_Covariates.Plot<-ggplot(QC.Dataframe,aes(x=nFeature_RNA, y=nCount_RNA ))+
geom_point(size=1)+
geom_hline(yintercept=minCount,color="blue",size=1,alpha=0.5)+
geom_hline(yintercept=maxCount,color="blue",size=1,alpha=0.5)+
geom_vline(xintercept=minFeatures,color="blue",size=1,alpha=0.5)+
geom_vline(xintercept=maxFeatures,color="blue",size=1,alpha=0.5)+
scale_x_continuous(trans='log2')+
scale_y_continuous(trans='log2')

ggsave2("Principal_QC_Covariates.Plot.png", Principal_QC_Covariates.Plot,units="in",width=6,height=3)





#### remove any nCount vs nFeature outliers. I think cells that deviate are probably techinical issues

nCount_nFeature_gam<-gam(nCount_RNA~poly(nFeature_RNA,2),data= QC.Dataframe)

QC.Dataframe$Count_Feature_gam_residual<-nCount_nFeature_gam$residuals

Count_Feature_gam_residual.StandardDeviation<-sd(QC.Dataframe$Count_Feature_gam_residual)
Count_Feature_gam_residual.Mean<-mean(QC.Dataframe$Count_Feature_gam_residual)
NumberOfSDs<-3



SelectingResidualCutoffs.Plot<-ggplot(QC.Dataframe,aes(x= Count_Feature_gam_residual))+
geom_histogram(bins=70)+
geom_vline(xintercept=Count_Feature_gam_residual.Mean)+
geom_vline(xintercept=Count_Feature_gam_residual.Mean+NumberOfSDs*Count_Feature_gam_residual.StandardDeviation,color="red")+
geom_vline(xintercept=Count_Feature_gam_residual.Mean-NumberOfSDs*Count_Feature_gam_residual.StandardDeviation,color="red")

ggsave2("SelectingResidualCutoffs.Plot.png", SelectingResidualCutoffs.Plot,units="in",width=6,height=3)

QC.Dataframe$nCount_nFeature_Deviation_Pass<-QC.Dataframe$Count_Feature_gam_residual > Count_Feature_gam_residual.Mean-NumberOfSDs*Count_Feature_gam_residual.StandardDeviation & QC.Dataframe$Count_Feature_gam_residual < Count_Feature_gam_residual.Mean+NumberOfSDs*Count_Feature_gam_residual.StandardDeviation

nFeature_nCount.Plot<-ggplot(QC.Dataframe,aes(x=nFeature_RNA,y=nCount_RNA))+
geom_point(aes(color=nCount_nFeature_Deviation_Pass))+
scale_color_manual(values=c("red","black"))+
geom_smooth(formula=y~poly(x,2),color="darkgrey")+
theme(legend.position = "none")

ggsave2("nFeature_nCount.Plot.png", nFeature_nCount.Plot,units="in",width=6,height=3)



#### remove cells beyond nCount and nFeature thresholds

nFeature_nCount.Cutoffs.Plot<-nFeature_nCount.Plot +
geom_hline(yintercept=minCount,color="blue",size=1,alpha=0.5)+
geom_hline(yintercept=maxCount,color="blue",size=1,alpha=0.5)+
geom_vline(xintercept=minFeatures,color="blue",size=1,alpha=0.5)+
geom_vline(xintercept=maxFeatures,color="blue",size=1,alpha=0.5)

ggsave2("nFeature_nCount.Cutoffs.Plot.png", nFeature_nCount.Cutoffs.Plot,units="in",width=6,height=3)


QC.Dataframe$nCount_nFeature_Threshold_Pass<-QC.Dataframe$nFeature_RNA > minFeatures & QC.Dataframe$nFeature_RNA < maxFeatures & QC.Dataframe$nCount_RNA > minCount & QC.Dataframe$nCount_RNA < maxCount


### No hard threshold on features
#QC.Dataframe.Passed<-subset(QC.Dataframe, nCount_nFeature_Deviation_Pass & fraction_mito_Pass)
QC.Dataframe.Passed<-subset(QC.Dataframe, nCount_nFeature_Threshold_Pass & nCount_nFeature_Deviation_Pass)



nrow(QC.Dataframe)

nrow(QC.Dataframe.Passed)


table(QC.Dataframe$nCount_nFeature_Threshold_Pass)

table(QC.Dataframe$nCount_nFeature_Deviation_Pass)





##### Filter based on QC

data.seurat<-subset(data.seurat, subset = orig.ident %in% QC.Dataframe.Passed$orig.ident)





QC.Panel<-plot_grid(nCount.Plot,nFeatureDistribution.Plot,SelectingResidualCutoffs.Plot,nFeature_nCount.Cutoffs.Plot,labels = "AUTO")

ggsave2("QC.Panel.png", QC.Panel,units="in",width=6,height=6)





data.seurat <- NormalizeData(data.seurat, normalization.method = "LogNormalize", scale.factor = 10000)



all.genes <- rownames(data.seurat)
data.seurat <- ScaleData(data.seurat, features = all.genes)



data.seurat <- FindVariableFeatures(object = data.seurat)
data.seurat <- RunPCA(data.seurat, features = NULL)


data.seurat <- FindNeighbors(data.seurat, dims = 1:10)



data.seurat <- FindClusters(data.seurat, resolution = 0.06)
data.seurat <- FindClusters(data.seurat, resolution = 0.07)
data.seurat <- FindClusters(data.seurat, resolution = 0.08)
data.seurat <- FindClusters(data.seurat, resolution = 0.09)
data.seurat <- FindClusters(data.seurat, resolution = 0.1)
data.seurat <- FindClusters(data.seurat, resolution = 0.11)
data.seurat <- FindClusters(data.seurat, resolution = 0.12)
data.seurat <- FindClusters(data.seurat, resolution = 0.13)
data.seurat <- FindClusters(data.seurat, resolution = 0.14)
data.seurat <- FindClusters(data.seurat, resolution = 0.15)
data.seurat <- FindClusters(data.seurat, resolution = 0.16)
data.seurat <- FindClusters(data.seurat, resolution = 0.17)
data.seurat <- FindClusters(data.seurat, resolution = 0.18)


#install.packages("clustree")
library(clustree)
ClusterPlot<-clustree(data.seurat)

ClusterPlot

ggsave("ClusterPlot.png",ClusterPlot,units="in",width = 12, height = 8)



data.seurat <- FindClusters(data.seurat, resolution = 0.15)




data.seurat <- RunUMAP(data.seurat, dims = 1:10)

UMAP<-DimPlot(data.seurat, reduction = "umap")

ggsave2("UMAP.png", UMAP,width=6,height=5,units="in")


ClusterPlot.Labeled<-ClusterPlot+
  geom_rect(aes(xmin = -7, xmax = 9, ymin = 2.5, ymax = 3.5), fill = "transparent", color = "black", size = 1,linetype="dashed")+
  scale_color_viridis(discrete = TRUE, option = "plasma", begin = 0.3, end = 0.9)


ClusterSelectionPanel<-plot_grid(ClusterPlot.Labeled,plot_grid(NULL,UMAP+ggtitle("RNA_snn_res = 0.15")+scale_color_viridis(discrete = TRUE, option = "plasma", begin = 0, end = 0.9),NULL,rel_heights = c(0.3,2,0.3),nrow=3),ncol=2,rel_widths = c(5,3.5), labels = "AUTO")

ggsave2("ClusterSelectionPanel.png", ClusterSelectionPanel,width=21.5,height=10,units="in")




markers <- FindAllMarkers(data.seurat, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, test.use = "MAST")
markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC)


top10 <- markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
TopGenes.Heatmap<-DoHeatmap(data.seurat, features = top10$gene) + NoLegend()

ggsave2("TopGenes.Heatmap.png", TopGenes.Heatmap,width=12,height=12,units="in")


write.table(subset(markers, cluster==0 & avg_logFC>1),file="Cluster.0.Markers.txt",quote=FALSE,sep='\t',row.names=FALSE)
write.table(subset(markers, cluster==1 & avg_logFC>1),file="Cluster.1.Markers.txt",quote=FALSE,sep='\t',row.names=FALSE)
write.table(subset(markers, cluster==2 & avg_logFC>1),file="Cluster.2.Markers.txt",quote=FALSE,sep='\t',row.names=FALSE)
write.table(subset(markers, cluster==3 & avg_logFC>1),file="Cluster.3.Markers.txt",quote=FALSE,sep='\t',row.names=FALSE)
write.table(subset(markers, cluster==4 & avg_logFC>1),file="Cluster.4.Markers.txt",quote=FALSE,sep='\t',row.names=FALSE)
write.table(subset(markers, cluster==5 & avg_logFC>1),file="Cluster.5.Markers.txt",quote=FALSE,sep='\t',row.names=FALSE)
write.table(subset(markers, cluster==6 & avg_logFC>1),file="Cluster.6.Markers.txt",quote=FALSE,sep='\t',row.names=FALSE)
write.table(subset(markers, cluster==7 & avg_logFC>1),file="Cluster.7.Markers.txt",quote=FALSE,sep='\t',row.names=FALSE)



FeaturePlot(data.seurat, features = c("Adipoq"))

DefiningCellTypes<-plot_grid(
  FeaturePlot(data.seurat, features = c("Adipoq"), order = TRUE)+theme(plot.title = element_text(face = "italic") ),
  FeaturePlot(data.seurat, features = c("Pdgfra"), order = TRUE)+theme(plot.title = element_text(face = "italic") ),
  FeaturePlot(data.seurat, features = c("Csf1r"), order = TRUE)+theme(plot.title = element_text(face = "italic") ),
  FeaturePlot(data.seurat, features = c("Cdh5"), order = TRUE)+theme(plot.title = element_text(face = "italic") ),
  FeaturePlot(data.seurat, features = c("Acta2"), order = TRUE)+theme(plot.title = element_text(face = "italic") ),
  FeaturePlot(data.seurat, features = c("Ptprc"), order = TRUE)+theme(plot.title = element_text(face = "italic") ),
  FeaturePlot(data.seurat, features = c("Cd2"), order = TRUE)+theme(plot.title = element_text(face = "italic") ),
  ncol = 4,
  labels = "AUTO"
)



new.cluster.ids <- c("MSC","Macrophage","VEC","VSMC","T cell","VEC","B cell","MSC.Diff")

names(new.cluster.ids) <- levels(data.seurat)


data.seurat <- RenameIdents(data.seurat, new.cluster.ids)
DimPlot(data.seurat, reduction = "umap", label = TRUE, pt.size = 1, label.size=4.5, repel=FALSE) + NoLegend()

new.cluster.id.order <- c("MSC","MSC.Diff","VEC","VSMC","Macrophage","T cell","B cell")

data.seurat@active.ident <- factor(x = data.seurat@active.ident, levels = new.cluster.id.order)



UMAP.Labeled<-DimPlot(data.seurat, reduction = "umap", label = TRUE, pt.size = 1, label.size=3, repel=FALSE) + NoLegend()
ggsave2("UMAP.Labeled.png", UMAP.Labeled,width=6,height=6,units="in")



DefiningCellTypes<-plot_grid(
  FeaturePlot(data.seurat, features = c("Adipoq"), order = TRUE)+theme(plot.title = element_text(face = "italic") ),
  FeaturePlot(data.seurat, features = c("Pdgfra"), order = TRUE)+theme(plot.title = element_text(face = "italic") ),
  FeaturePlot(data.seurat, features = c("Csf1r"), order = TRUE)+theme(plot.title = element_text(face = "italic") ),
  FeaturePlot(data.seurat, features = c("Cdh5"), order = TRUE)+theme(plot.title = element_text(face = "italic") ),
  FeaturePlot(data.seurat, features = c("Acta2"), order = TRUE)+theme(plot.title = element_text(face = "italic") ),
  FeaturePlot(data.seurat, features = c("Ptprc"), order = TRUE)+theme(plot.title = element_text(face = "italic") ),
  FeaturePlot(data.seurat, features = c("Cd2"), order = TRUE)+theme(plot.title = element_text(face = "italic") ),
  UMAP.Labeled+ggtitle("Cell Type")+theme(plot.title = element_text(face = "italic") ),
  ncol = 4,
  labels = "AUTO"
)
ggsave2("DefiningCellTypes.png",DefiningCellTypes,width=12,height=6)




ListOfGenes<-c("F2r","Plcd1","Hexb","Nnat","Cdkn1c","Hmgcr","Car3","Tnfrsf11a","Srgn","Adipoq","Pdgfra")




Data.DE<-FindMarkers(data.seurat, ident.1 = "MSC", ident.2 = "MSC.Diff", min.pct = 0, test.use = "MAST", logfc.threshold=0)


Data.DE$Gene<-rownames(Data.DE)
rownames(Data.DE)<-NULL

Data.DE$GOI<-Data.DE$Gene %in% ListOfGenes

Data.DE$pct.change<-Data.DE$pct.1-Data.DE$pct.2
Data.DE$pct.change.scaled<- unlist(lapply(1:nrow(Data.DE),function(ROW) (Data.DE$pct.1[ROW]-Data.DE$pct.2[ROW])/max(c(Data.DE$pct.1[ROW],Data.DE$pct.2[ROW]))))












Contributions_InVitro<-read.delim("Contributions.txt",sep="\t", header = TRUE)
Contributions_InVivo<-read.delim("Contributions.Literature.txt",sep="\t", header = TRUE)

Contributions_InVitro$CellType<-factor(Contributions_InVitro$CellType,levels = c("MSC","MSC.Diff","VSMC","VEC","B cell","Macrophage","T cell"))
Contributions_InVivo$CellType<-factor(Contributions_InVivo$CellType,levels = c("MSC","MSC.Diff","VSMC","VEC","B cell","Macrophage","T cell"))

Data.DE$AdipocyteContribution.InVitro<-unlist(lapply(Data.DE$Gene,function(GENE) sum(subset(Contributions_InVitro,Gene==GENE & CellType %in% c("MSC","MSC.Diff"))$FractionOfTotal,na.rm = TRUE)))
Data.DE$AdipocyteContribution.InVivo<-unlist(lapply(Data.DE$Gene,function(GENE) sum(subset(Contributions_InVivo,Gene==GENE & CellType %in% c("MSC","MSC.Diff"))$FractionOfTotal,na.rm = TRUE)))

InVivoThreshold<-0.7
InVitroThreshold<-0.3

Data.DE$InVivoTest<-Data.DE$AdipocyteContribution.InVivo>=InVivoThreshold
Data.DE$InVitroTest<-Data.DE$AdipocyteContribution.InVitro>=InVitroThreshold



pct.change.scaled.threshold<-0.4
avg.logFC.threshold<-0.3


VolcanoPlot<-ggplot(Data.DE %>% arrange(GOI) ,aes(x=avg_logFC,y=-log10(p_val_adj),color=GOI))+
  geom_point(size=0.5)+
  theme(legend.position = "none")+
  scale_color_manual(values = c("grey","blue"))+
  geom_vline(xintercept = avg.logFC.threshold)+
  geom_vline(xintercept = -avg.logFC.threshold)+
  geom_hline(yintercept = -log10(0.0001))+
  ylab("-Log10(Bonferoni)")+
  xlab("Mean Log2 Fold-Change")

ggsave(filename = "VolcanoPlot.png",plot = VolcanoPlot, units = "in", width = 3, height = 3)


PercentShiftFoldChangePlot<-ggplot(Data.DE %>% arrange(GOI) ,aes(x=pct.change.scaled,y=avg_logFC,color=GOI))+
  geom_point(size=0.5)+
  theme(legend.position = "none")+
  scale_color_manual(values = c("grey","blue"))+
  geom_hline(yintercept = avg.logFC.threshold)+
  geom_hline(yintercept = -avg.logFC.threshold)+
  geom_vline(xintercept = pct.change.scaled.threshold)+
  geom_vline(xintercept = -pct.change.scaled.threshold)+
  xlab("Scaled Change in Percent Expressing")+
  ylab("Mean Log2 Fold-Change")

ggsave(filename = "PercentShiftFoldChangePlot.png",plot = PercentShiftFoldChangePlot, units = "in", width = 3, height = 3)





DE.Panel<-plot_grid(UMAP.Labeled,VolcanoPlot,PercentShiftFoldChangePlot,ncol=3, labels = "AUTO")
ggsave(filename = "DE.Panel.png",plot = DE.Panel,units = "in",width = 9, height = 3 )




Data.DE$pct.change.scaled.Test<-abs(Data.DE$pct.change.scaled)>=pct.change.scaled.threshold
Data.DE$avg_logFC.Test<-abs(Data.DE$avg_logFC)>=avg.logFC.threshold


##### With using cell type specific thresholds
nrow(subset(Data.DE, InVivoTest & InVitroTest & ( pct.change.scaled.Test | avg_logFC.Test ) ))
### 6717

nrow(subset(Data.DE, InVivoTest & InVitroTest & ( pct.change.scaled.Test | avg_logFC.Test ) ))/nrow(Data.DE)
### 0.5687071



##### Without using cell type specific thresholds
nrow(subset(Data.DE,( pct.change.scaled.Test | avg_logFC.Test ) ))
### 9215

nrow(subset(Data.DE, ( pct.change.scaled.Test | avg_logFC.Test ) ))/nrow(Data.DE)
### 0.7802049




write.table(x = Data.DE,file = "All.Adipocyre.Data.DE.txt", sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)

subset(Data.DE, GOI)







ASC.data<-subset(data.seurat,cells=names(Idents(data.seurat))[Idents(data.seurat) %in% c("MSC","MSC.Diff") ])


all.genes <- rownames(ASC.data)
ASC.data <- ScaleData(ASC.data, features = all.genes)



ASC.data <- FindVariableFeatures(object = ASC.data)
ASC.data <- RunPCA(ASC.data, features = NULL)


ASC.data <- FindNeighbors(ASC.data, dims = 1:10)




ASC.data <- FindClusters(ASC.data, resolution = 0.01)
ASC.data <- FindClusters(ASC.data, resolution = 0.02)
ASC.data <- FindClusters(ASC.data, resolution = 0.03)
ASC.data <- FindClusters(ASC.data, resolution = 0.04)
ASC.data <- FindClusters(ASC.data, resolution = 0.05)
ASC.data <- FindClusters(ASC.data, resolution = 0.06)
ASC.data <- FindClusters(ASC.data, resolution = 0.07)
ASC.data <- FindClusters(ASC.data, resolution = 0.08)
ASC.data <- FindClusters(ASC.data, resolution = 0.09)
ASC.data <- FindClusters(ASC.data, resolution = 0.1)
ASC.data <- FindClusters(ASC.data, resolution = 0.11)
ASC.data <- FindClusters(ASC.data, resolution = 0.12)
ASC.data <- FindClusters(ASC.data, resolution = 0.13)
ASC.data <- FindClusters(ASC.data, resolution = 0.14)
ASC.data <- FindClusters(ASC.data, resolution = 0.15)
ASC.data <- FindClusters(ASC.data, resolution = 0.16)
ASC.data <- FindClusters(ASC.data, resolution = 0.17)
ASC.data <- FindClusters(ASC.data, resolution = 0.18)
ASC.data <- FindClusters(ASC.data, resolution = 0.19)
ASC.data <- FindClusters(ASC.data, resolution = 0.20)
ASC.data <- FindClusters(ASC.data, resolution = 0.21)
ASC.data <- FindClusters(ASC.data, resolution = 0.22)
ASC.data <- FindClusters(ASC.data, resolution = 0.23)
ASC.data <- FindClusters(ASC.data, resolution = 0.24)
ASC.data <- FindClusters(ASC.data, resolution = 0.25)


ClusterPlot.ASC<-clustree(ASC.data)



ggsave("ClusterPlot.ASC.png",ClusterPlot.ASC,units="in",width = 12, height = 10)




#### Find cluster that minimizes Adipoq var

RES.NAME.LIST<-paste("RNA_snn_res",seq(0.01,0.25,0.01),sep=".")

minimize.SD.DF<-data.frame(
  Resolution=seq(0.01,0.25,0.01),
  mean.sd=unlist(lapply(
  RES.NAME.LIST,
  function(RES.NAME) mean(unlist(lapply(1:length(unique(ASC.data@meta.data[[RES.NAME]])), function(CLUSTERINDEX) sd(ASC.data@assays$RNA@data[rownames(ASC.data@assays$RNA@data)=="Adipoq",ASC.data@meta.data[[RES.NAME]]==CLUSTERINDEX-1][ASC.data@assays$RNA@data[rownames(ASC.data@assays$RNA@data)=="Adipoq",ASC.data@meta.data[[RES.NAME]]==CLUSTERINDEX-1]>0]) )))
    )),
  mean.exp.sd=unlist(lapply(
    RES.NAME.LIST,
    function(RES.NAME) mean(unlist(lapply(1:length(unique(ASC.data@meta.data[[RES.NAME]])), function(CLUSTERINDEX) sd(ASC.data@assays$RNA@data[rownames(ASC.data@assays$RNA@data)=="Adipoq",ASC.data@meta.data[[RES.NAME]]==CLUSTERINDEX-1]>0) )))
  )),
  NumberOfClusters=unlist(lapply(
    RES.NAME.LIST,
    function(RES.NAME) length(unique(ASC.data@meta.data[[RES.NAME]]))
  ))
)


write.table(minimize.SD.DF, file = "minimize.SD.DF.txt",sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)


ASC.data <- FindClusters(ASC.data, resolution = 0.17)

ASC.data <- RunUMAP(ASC.data, dims = 1:10)

ASC.UMAP<-DimPlot(ASC.data, reduction = "umap", label = TRUE, pt.size = 1) + NoLegend()


ggsave2("ASC.UMAP.png", ASC.UMAP,width=6,height=5,units="in")




AdipoqClusterPanel.A<-FeaturePlot(ASC.data,features = "Adipoq", order = TRUE, pt.size = 1)

AdipoqClusterPanel.B<-ASC.UMAP

AdipoqClusterPanel.C<-ggplot(minimize.SD.DF,aes(x=Resolution,y=mean.sd))+
  geom_vline(xintercept = 0.17, color="red", size=1)+
  geom_line()+
  ylab("Mean Count Variance")

AdipoqClusterPanel.D<-ggplot(minimize.SD.DF,aes(x=Resolution,y=mean.exp.sd))+
  geom_vline(xintercept = 0.17, color="red", size=1)+
  geom_line()+
  ylab("Percent Expressing Variance")

AdipoqClusterPanel.Left<-plot_grid(AdipoqClusterPanel.A,AdipoqClusterPanel.B,AdipoqClusterPanel.C,AdipoqClusterPanel.D,nrow=4,labels = c("A","B","C","D"))

AdipoqClusterPanel<-plot_grid(AdipoqClusterPanel.Left,ClusterPlot.ASC, nrow=1, rel_widths = c(1,1.5), rel_heights = c(1,1), labels = c("","E"))

ggsave(filename = "AdipoqClusterPanel.png", plot = AdipoqClusterPanel, units = "in", width = 8, height = 12)








new.cluster.ids.2 <- c("1","0","0","2","X")
names(new.cluster.ids.2) <- levels(ASC.data)
ASC.data <- RenameIdents(ASC.data, new.cluster.ids.2)


UMAP.Labeled.ASC<-DimPlot(ASC.data, reduction = "umap", label = TRUE, pt.size = 2) + NoLegend()

ggsave2("UMAP.Labeled.ASC.png", UMAP.Labeled.ASC,width=6,height=6,units="in")


ASC.Pre<-subset(ASC.data,cells=names(Idents(ASC.data))[Idents(ASC.data) %in% c(0) ])
ASC.Mid<-subset(ASC.data,cells=names(Idents(ASC.data))[Idents(ASC.data) %in% c(1) ])
ASC.Diff<-subset(ASC.data,cells=names(Idents(ASC.data))[Idents(ASC.data) %in% c(2) ])


UMAP.Panel<-plot_grid(
FeaturePlot(ASC.data, features = c("F2r"), order=TRUE)+theme(axis.text=element_text(size=5), axis.title=element_text(size=9), legend.text = element_text(size = 6), legend.key.size = unit(0.25, 'cm'))+theme(plot.title = element_text(face = "italic") ),
FeaturePlot(ASC.data, features = c("Plcd1"), order=TRUE)+theme(axis.text=element_text(size=5), axis.title=element_text(size=9), legend.text = element_text(size = 6), legend.key.size = unit(0.25, 'cm'))+theme(plot.title = element_text(face = "italic") ),
FeaturePlot(ASC.data, features = c("Hexb"), order=TRUE)+theme(axis.text=element_text(size=5), axis.title=element_text(size=9), legend.text = element_text(size = 6), legend.key.size = unit(0.25, 'cm'))+theme(plot.title = element_text(face = "italic") ),
FeaturePlot(ASC.data, features = c("Nnat"), order=TRUE)+theme(axis.text=element_text(size=5), axis.title=element_text(size=9), legend.text = element_text(size = 6), legend.key.size = unit(0.25, 'cm'))+theme(plot.title = element_text(face = "italic") ),
FeaturePlot(ASC.data, features = c("Cdkn1c"), order=TRUE)+theme(axis.text=element_text(size=5), axis.title=element_text(size=9), legend.text = element_text(size = 6), legend.key.size = unit(0.25, 'cm'))+theme(plot.title = element_text(face = "italic") ),
FeaturePlot(ASC.data, features = c("Hmgcr"), order=TRUE)+theme(axis.text=element_text(size=5), axis.title=element_text(size=9), legend.text = element_text(size = 6), legend.key.size = unit(0.25, 'cm'))+theme(plot.title = element_text(face = "italic") ),
FeaturePlot(ASC.data, features = c("Car3"), order=TRUE)+theme(axis.text=element_text(size=5), axis.title=element_text(size=9), legend.text = element_text(size = 6), legend.key.size = unit(0.25, 'cm'))+theme(plot.title = element_text(face = "italic") ),
FeaturePlot(ASC.data, features = c("Srgn"), order=TRUE)+theme(axis.text=element_text(size=5), axis.title=element_text(size=9), legend.text = element_text(size = 6), legend.key.size = unit(0.25, 'cm'))+theme(plot.title = element_text(face = "italic") ),
FeaturePlot(ASC.data, features = c("Adipoq"), order=TRUE)+theme(axis.text=element_text(size=5), axis.title=element_text(size=9), legend.text = element_text(size = 6), legend.key.size = unit(0.25, 'cm'))+theme(plot.title = element_text(face = "italic") ),
UMAP.Labeled.ASC+theme(plot.title = element_text(face = "italic") ),
labels = "AUTO",
ncol = 5
)


ggsave("UMAP.Panel.png", UMAP.Panel, width=12, height=4, units='in')




BLENDCOLS<-c("grey93",viridis(1,begin=0.5,option="viridis"),viridis(1,begin=0.2,option="viridis"))
Rev.BLENDCOLS<-c("grey93",viridis(1,begin=0.2,option="magma"),viridis(1,begin=0.5,option="magma"))



Adipoq.Pdgfra.Split<-FeaturePlot(ASC.data, features = c("Adipoq","Pdgfra"), order=TRUE, pt.size = 0.5, blend=TRUE,cols=Rev.BLENDCOLS, combine = FALSE)
F2r.Nnat.Split<-FeaturePlot(ASC.data, features = c("F2r","Nnat"), order=TRUE, pt.size = 0.5, blend=TRUE,cols=BLENDCOLS, combine = FALSE)




TrajectoryPanel<-plot_grid(
  plot_grid(
  Adipoq.Pdgfra.Split[[1]]+NoLegend()+theme(plot.title = element_text(face = "italic") )+theme(axis.title=element_text(size=9)),
  Adipoq.Pdgfra.Split[[2]]+NoLegend()+theme(plot.title = element_text(face = "italic") )+theme(axis.title=element_text(size=9)),
  Adipoq.Pdgfra.Split[[3]]+NoLegend()+ggtitle("Adipoq:Pdgfra")+theme(plot.title = element_text(face = "italic") )+theme(axis.title=element_text(size=9)),
  F2r.Nnat.Split[[2]]+NoLegend()+theme(plot.title = element_text(face = "italic") )+theme(axis.title=element_text(size=9)),
  F2r.Nnat.Split[[1]]+NoLegend()+theme(plot.title = element_text(face = "italic") )+theme(axis.title=element_text(size=9)),
  F2r.Nnat.Split[[3]]+NoLegend()+ggtitle("F2r:Nnat")+theme(plot.title = element_text(face = "italic") )+theme(axis.title=element_text(size=9)),
  ncol=3,
  labels = c("A","B","C","E","F","G")
  ),
  plot_grid(
    plot_grid(
      Adipoq.Pdgfra.Split[[4]]+ggtitle("")+theme(plot.margin = margin(0,0,0,0))+theme(axis.text=element_text(size=7),axis.title=element_text(size=9)),
      F2r.Nnat.Split[[4]]+ggtitle("")+theme(plot.margin = margin(0,0,0,0))+theme(axis.text=element_text(size=7),axis.title=element_text(size=9)),
      labels=c("D","H"),
      ncol=1
    ),
    UMAP.Labeled.ASC+scale_color_viridis(option = "inferno", discrete = TRUE, begin = 1, end = 0.4)+theme(axis.title=element_text(size=9),plot.title = element_text(face = "italic"))+ggtitle("Clusters"),
    ncol=1,
    labels=c("","I"),
    rel_heights = c(1,1)
    ),
  rel_widths = c(3,1)
)


ggsave("TrajectoryPanel.png", TrajectoryPanel, width=9, height=4, units='in')



Adipoq.Pdgfra<-FeaturePlot(ASC.data, features = c("Adipoq","Pdgfra"), order=TRUE, pt.size = 0.5, blend=TRUE,cols=Rev.BLENDCOLS)+theme(axis.text=element_text(size=5), axis.title=element_text(size=9), legend.text = element_text(size = 6), legend.key.size = unit(0.25, 'cm'),plot.title = element_blank())+NoLegend()


ggsave("Adipoq.Pdgfra.png", Adipoq.Pdgfra, width=9, height=2.5, units='in')


# F0/1, F1 Net, GTEX Visc, GTEX Sub

# Neg, Neg, Neg, Neg
F2r.Nnat<-FeaturePlot(ASC.data, features = c("F2r","Nnat"), order=TRUE, pt.size = 0.5, blend=TRUE,cols=BLENDCOLS)+theme(axis.text=element_text(size=5), axis.title=element_text(size=9), legend.text = element_text(size = 6), legend.key.size = unit(0.25, 'cm'),plot.title = element_blank())+NoLegend()
# Neg

ggsave("F2r.Nnat.png", F2r.Nnat, width=9, height=2.5, units='in')

# Neg, Neg, Neg, Neg
Hexb.Nnat<-FeaturePlot(ASC.data, features = c("Hexb","Nnat"), order=TRUE, pt.size = 0.5, blend=TRUE,cols=BLENDCOLS)+theme(axis.text=element_text(size=5), axis.title=element_text(size=9), legend.text = element_text(size = 6), legend.key.size = unit(0.25, 'cm'),plot.title = element_blank())+NoLegend()
# Neg

ggsave("Hexb.Nnat.png", Hexb.Nnat, width=9, height=2.5, units='in')


# Neg, Neg, Pos, No
Hmgcr.Nnat<-FeaturePlot(ASC.data, features = c("Hmgcr","Nnat"), order=TRUE, pt.size = 0.5, blend=TRUE,cols=BLENDCOLS)+theme(axis.text=element_text(size=5), axis.title=element_text(size=9), legend.text = element_text(size = 6), legend.key.size = unit(0.25, 'cm'),plot.title = element_blank())+NoLegend()
# Neg

ggsave("Hmgcr.Nnat.png", Hmgcr.Nnat, width=9, height=2.5, units='in')



# No, Neg, Neg, Neg
F2r.Cdkn1c<-FeaturePlot(ASC.data, features = c("F2r","Cdkn1c"), order=TRUE, pt.size = 0.5, blend=TRUE,cols=BLENDCOLS)+theme(axis.text=element_text(size=5), axis.title=element_text(size=9), legend.text = element_text(size = 6), legend.key.size = unit(0.25, 'cm'),plot.title = element_blank())+NoLegend()
# Pos

ggsave("F2r.Cdkn1c.png", F2r.Cdkn1c, width=9, height=2.5, units='in')


# No, Neg, Pos, Neg
Hmgcr.Cdkn1c<-FeaturePlot(ASC.data, features = c("Hmgcr","Cdkn1c"), order=TRUE, pt.size = 0.5, blend=TRUE,cols=BLENDCOLS)+theme(axis.text=element_text(size=5), axis.title=element_text(size=9), legend.text = element_text(size = 6), legend.key.size = unit(0.25, 'cm'),plot.title = element_blank())+NoLegend()
# Pos

ggsave("Hmgcr.Cdkn1c.png", Hmgcr.Cdkn1c, width=9, height=2.5, units='in')

# Neg, Neg, Pos, No
Srgn.Cdkn1c<-FeaturePlot(ASC.data, features = c("Srgn","Cdkn1c"), order=TRUE, pt.size = 0.5, blend=TRUE,cols=BLENDCOLS)+theme(axis.text=element_text(size=5), axis.title=element_text(size=9), legend.text = element_text(size = 6), legend.key.size = unit(0.25, 'cm'),plot.title = element_blank())+NoLegend()
# Pos

ggsave("Srgn.Cdkn1c.png", Srgn.Cdkn1c, width=9, height=2.5, units='in')



# Pos, Pos, Pos, No
Car3.Plcd1<-FeaturePlot(ASC.data, features = c("Car3","Plcd1"), order=TRUE, pt.size = 0.5, blend=TRUE,cols=BLENDCOLS)+theme(axis.text=element_text(size=5), axis.title=element_text(size=9), legend.text = element_text(size = 6), legend.key.size = unit(0.25, 'cm'),plot.title = element_blank())+NoLegend()
# Neg

ggsave("Car3.Plcd1.png", Car3.Plcd1, width=9, height=2.5, units='in')








ASC.data@active.ident<-factor(ASC.data@active.ident,levels = c(0,1,2,"X"))



UMAP.ASC<-DimPlot(ASC.data, reduction = "umap", label = TRUE, pt.size = 0.5, label.size=4) + NoLegend() + theme(axis.text=element_text(size=5), axis.title=element_text(size=9), legend.text = element_text(size = 6), legend.key.size = unit(0.25, 'cm'))


ggsave("UMAP.ASC.png", UMAP.ASC, width=2, height=2, units='in')

AdipoqUMAP<-FeaturePlot(ASC.data, features = c("Adipoq"), order=TRUE, pt.size = 0.5)+theme(axis.text=element_text(size=5), axis.title=element_text(size=9), legend.text = element_text(size = 6), legend.key.size = unit(0.25, 'cm'))+NoLegend()

blankPlot <- ggplot()+geom_blank(aes(1,1)) + cowplot::theme_nothing()

Panel.DotPlot<-DotPlot(ASC.data, features = ListOfGenes, scale=TRUE) + RotatedAxis() + ylab("") + theme( legend.text = element_text(size = 6), legend.key.size = unit(0.25, 'cm'), legend.title = element_text(size = 6))

DotplotLegend<-get_legend(Panel.DotPlot)

Panel.DotPlot<-Panel.DotPlot+theme(legend.position = "none")

GROB.Leg<-arrangeGrob(blankPlot, blankPlot, blankPlot, DotplotLegend, blankPlot ,nrow=3, heights=c(2,0.1,2), widths=c(1,30), ncol=2)

grid.arrange(GROB.Leg)


GROB<-arrangeGrob(Panel.DotPlot, GROB.Leg, ncol=2, widths=c(5,1))

grid.arrange(GROB)


GROB2<-arrangeGrob(UMAP.ASC, AdipoqUMAP,ncol=2,widths=c(2,2))

DotPlotPanel<-grid.arrange(GROB2, GROB, nrow=2)



ggsave("DotPlotPanel.png", DotPlotPanel, width=5, height=5, units='in')



Vln.F2r<-VlnPlot(ASC.data, features = ListOfGenes[1], pt.size=0.1)+ylim(0.0001,2.5)+theme(legend.position = "none", plot.title = element_text(face = "italic") )+xlab("Subtype")
Vln.Plcd1<-VlnPlot(ASC.data, features = ListOfGenes[2], pt.size=0.1)+ylim(0.0001,3)+theme(legend.position = "none", plot.title = element_text(face = "italic") )+xlab("Subtype")
Vln.Hexb<-VlnPlot(ASC.data, features = ListOfGenes[3], pt.size=0.1)+ylim(0.0001,3.5)+theme(legend.position = "none", plot.title = element_text(face = "italic") )+xlab("Subtype")
Vln.Nnat<-VlnPlot(ASC.data, features = ListOfGenes[4], pt.size=0.1)+ylim(0.0001,3.5)+theme(legend.position = "none", plot.title = element_text(face = "italic") )+xlab("Subtype")
Vln.Cdkn1c<-VlnPlot(ASC.data, features = ListOfGenes[5], pt.size=0.1)+ylim(0.0001,3.5)+theme(legend.position = "none", plot.title = element_text(face = "italic") )+xlab("Subtype")
Vln.Hmgcr<-VlnPlot(ASC.data, features = ListOfGenes[6], pt.size=0.1)+ylim(0.0001,2.5)+theme(legend.position = "none", plot.title = element_text(face = "italic") )+xlab("Subtype")
Vln.Car3<-VlnPlot(ASC.data, features = ListOfGenes[7], pt.size=0.1)+ylim(0.0001,6.5)+theme(legend.position = "none", plot.title = element_text(face = "italic") )+xlab("Subtype")
Vln.Tnfrsf11a<-VlnPlot(ASC.data, features = ListOfGenes[8], pt.size=0.1)+ylim(0.0001,2)+theme(legend.position = "none", plot.title = element_text(face = "italic") )+xlab("Subtype")
Vln.Srgn<-VlnPlot(ASC.data, features = ListOfGenes[9], pt.size=0.1)+ylim(0.0001,3.5)+theme(legend.position = "none", plot.title = element_text(face = "italic") )+xlab("Subtype")
Vln.Adipoq<-VlnPlot(ASC.data, features = ListOfGenes[10], pt.size=0.1)+ylim(0.0001,5.5)+theme(legend.position = "none", plot.title = element_text(face = "italic") )+xlab("Subtype")
Vln.Pdgfra<-VlnPlot(ASC.data, features = ListOfGenes[11], pt.size=0.1)+ylim(0.0001,4)+theme(legend.position = "none", plot.title = element_text(face = "italic") )+xlab("Subtype")


Vln.A<-plot_grid(Vln.Pdgfra,Vln.F2r,Vln.Plcd1,Vln.Hexb,Vln.Cdkn1c,Vln.Hmgcr, nrow=1, labels = c("G","H","I","J","K","L"))

Vln.B<-plot_grid(Vln.Adipoq,Vln.Nnat,Vln.Car3,nrow=1, labels = c("B","C","D"))

Vln.C<-plot_grid(Vln.Tnfrsf11a,Vln.Srgn,nrow=1, labels = c("E","F"))


Proadip.Title<-ggplot()+
  annotate("text",x=1,y=1,size=8,label='bold("Upregulated with differentiation")', parse=TRUE)+
  theme_void()

Vln.ProAdip<-plot_grid(Proadip.Title,Vln.B,rel_heights = c(0.2,1), nrow=2)


Immune.Title<-ggplot()+
  annotate("text",x=1,y=1,size=8,label='bold("Immune")', parse=TRUE)+
  theme_void()

Vln.Immune<-plot_grid(Immune.Title,Vln.C,rel_heights = c(0.2,1), nrow=2)


AntiAdip.Title<-ggplot()+
  annotate("text",x=1,y=1,size=8,label='bold("Downregulated with differentiation")', parse=TRUE)+
  theme_void()

Vln.AntiAdip<-plot_grid(AntiAdip.Title,Vln.A,rel_heights = c(0.2,1), nrow=2)


UMAP.Title<-ggplot()+
  annotate("text",x=1,y=1,size=8,label='bold("Cell Types")', parse=TRUE)+
  theme_void()

Vln.UMAP<-plot_grid(UMAP.Title,plot_grid(UMAP.ASC,labels=c("A")),rel_heights = c(0.2,1), nrow=2)


Vln.Panel<-plot_grid(plot_grid(Vln.UMAP,NULL,Vln.ProAdip,NULL,Vln.Immune,nrow=1, rel_widths = c(1,0.2,3,0.2,2) ),plot_grid(NULL,Vln.AntiAdip,NULL,nrow = 1, rel_widths = c(0.2,6,0.2)),nrow=2)



ggsave(filename = "Vln.Panel.png",plot = Vln.Panel,units = "in",width = 16, height = 16*(5/12) )






Data.DE <- FindMarkers(ASC.data, ident.1 = "1", ident.2 = "2", min.cells.group = 1, min.cells.feature = 1, min.pct = 0, logfc.threshold = 0, only.pos = FALSE, test.use ="MAST")

Data.DE$Gene<-rownames(Data.DE)
rownames(Data.DE)<-NULL

Data.DE$GOI<-Data.DE$Gene %in% ListOfGenes

Data.DE$pct.change<-Data.DE$pct.1-Data.DE$pct.2
Data.DE$pct.change.scaled<- unlist(lapply(1:nrow(Data.DE),function(ROW) (Data.DE$pct.1[ROW]-Data.DE$pct.2[ROW])/max(c(Data.DE$pct.1[ROW],Data.DE$pct.2[ROW]))))





Data.DE$AdipocyteContribution.InVitro<-unlist(lapply(Data.DE$Gene,function(GENE) sum(subset(Contributions_InVitro,Gene==GENE & CellType %in% c("MSC","MSC.Diff"))$FractionOfTotal,na.rm = TRUE)))
Data.DE$AdipocyteContribution.InVivo<-unlist(lapply(Data.DE$Gene,function(GENE) sum(subset(Contributions_InVivo,Gene==GENE & CellType %in% c("MSC","MSC.Diff"))$FractionOfTotal,na.rm = TRUE)))

InVivoThreshold<-0.7
InVitroThreshold<-0.3

Data.DE$InVivoTest<-Data.DE$AdipocyteContribution.InVivo>=InVivoThreshold
Data.DE$InVitroTest<-Data.DE$AdipocyteContribution.InVitro>=InVitroThreshold



pct.change.scaled.threshold<-0.4
avg.logFC.threshold<-0.5


VolcanoPlot.ASC<-ggplot(Data.DE %>% arrange(GOI) ,aes(x=avg_logFC,y=-log10(p_val_adj),color=GOI))+
  geom_point(size=0.5)+
  theme(legend.position = "none")+
  scale_color_manual(values = c("grey","blue"))+
  geom_vline(xintercept = avg.logFC.threshold)+
  geom_vline(xintercept = -avg.logFC.threshold)+
  geom_hline(yintercept = -log10(0.0001))+
  ylab("-Log10(Bonferoni)")+
  xlab("Mean Log2 Fold-Change")

ggsave(filename = "VolcanoPlot.ASC.png",plot = VolcanoPlot.ASC, units = "in", width = 3, height = 3)


PercentShiftFoldChangePlot.ASC<-ggplot(Data.DE %>% arrange(GOI) ,aes(x=pct.change.scaled,y=avg_logFC,color=GOI))+
  geom_point(size=0.5)+
  theme(legend.position = "none")+
  scale_color_manual(values = c("grey","blue"))+
  geom_hline(yintercept = avg.logFC.threshold)+
  geom_hline(yintercept = -avg.logFC.threshold)+
  geom_vline(xintercept = pct.change.scaled.threshold)+
  geom_vline(xintercept = -pct.change.scaled.threshold)+
  xlab("Scaled Change in Percent Expressing")+
  ylab("Mean Log2 Fold-Change")

ggsave(filename = "PercentShiftFoldChangePlot.ASC.png",plot = PercentShiftFoldChangePlot, units = "in", width = 3, height = 3)


Data.DE$pct.change.scaled.Test<-abs(Data.DE$pct.change.scaled)>=pct.change.scaled.threshold
Data.DE$avg_logFC.Test<-abs(Data.DE$avg_logFC)>=avg.logFC.threshold


##### With using cell type specific thresholds
nrow(subset(Data.DE, InVivoTest & InVitroTest & ( pct.change.scaled.Test | avg_logFC.Test ) ))
### 6714

nrow(subset(Data.DE, InVivoTest & InVitroTest & ( pct.change.scaled.Test | avg_logFC.Test ) ))/nrow(Data.DE)
### 0.5693208



##### Without using cell type specific thresholds
nrow(subset(Data.DE,( pct.change.scaled.Test | avg_logFC.Test ) ))
### 9265

nrow(subset(Data.DE, ( pct.change.scaled.Test | avg_logFC.Test ) ))/nrow(Data.DE)
### 0.7856355




write.table(x = Data.DE,file = "ASC.Type1.Type2.Data.DE.txt", sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)




subset(Data.DE,GOI)



ASC.DE.Panel<-plot_grid(UMAP.Labeled.ASC,VolcanoPlot.ASC,PercentShiftFoldChangePlot.ASC,ncol=3, labels = "AUTO")
ggsave(filename = "ASC.DE.Panel.png",plot = ASC.DE.Panel,units = "in",width = 9, height = 3 )






saveRDS(data.seurat, file = "data.seurat")

saveRDS(ASC.data, file = "ASC.data")




