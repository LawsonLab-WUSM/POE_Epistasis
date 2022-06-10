
rm(list=ls())
library(dplyr)
library(patchwork)
library(cowplot)
library(ggplot2)
library(gridExtra)
library(viridis)
theme_set(theme_cowplot())

theme_set(theme(
  panel.background = element_rect(fill="white"),
  panel.grid.major=element_line(color="gray90", linetype="solid"), #adds major grid lines
  panel.grid.minor=element_line(color="gray95", linetype="solid", size=0.65), #adds minor grid lines	
  panel.border=element_rect(fill=NA, color="black", size=1.5, linetype="solid"), #draws a black border around plot
  axis.line=element_blank()
))







CellTypeBreakdown<-read.delim("CellTypeCompositions.txt",header=TRUE,sep="\t")
CellTypeBreakdown$CellType<-factor(CellTypeBreakdown$CellType,levels = c("MSC","MSC.Diff","VSMC","VEC","B cell","Macrophage","T cell"))

CellTypeBreakdown


CellTyepColors<-c(
  viridis(2,option = "plasma", begin = 0.075,end = 0.2),
  viridis(2,option = "plasma", begin = 0.45,end = 0.6),
  viridis(3,option = "plasma", begin = 0.8,end = 0.95)
)


Pie.InVitro<-ggplot(subset(CellTypeBreakdown, Source=="In vitro"), aes(x="",y=CellFraction,fill=CellType))+
  geom_bar(stat = "identity", width = 1, color="white")+
  coord_polar("y",start=0)+
  theme_void()+
  labs(fill="Cell type")+
  theme(legend.position = "none")+
  scale_fill_manual(values = CellTyepColors)+
  ggtitle(expression(italic("In Vitro")))




Pie.InVivo<-ggplot(subset(CellTypeBreakdown, Source=="In vivo"), aes(x="",y=CellFraction,fill=CellType))+
  geom_bar(stat = "identity", width = 1, color="white")+
  coord_polar("y",start=0)+
  theme_void()+
  labs(fill="Cell type")+
  theme(legend.position = "bottom",legend.box.margin = margin(5,0,5,0))+
  scale_fill_manual(values = CellTyepColors)+
  ggtitle(expression(italic("In Vivo")))

Pie.Legend<-get_legend(Pie.InVivo)

Pie.InVivo<-Pie.InVivo+theme(legend.position = "none")





Contributions_InVitro<-read.delim("Contributions.txt",sep="\t", header = TRUE)

Contributions_InVivo<-read.delim("Contributions.Literature.txt",sep="\t", header = TRUE)



Contributions_InVitro$CellType<-factor(Contributions_InVitro$CellType,levels = c("MSC","MSC.Diff","VSMC","VEC","B cell","Macrophage","T cell"))
Contributions_InVivo$CellType<-factor(Contributions_InVivo$CellType,levels = c("MSC","MSC.Diff","VSMC","VEC","B cell","Macrophage","T cell"))






sort(intersect(
  subset(Contributions_InVitro,CellType=="MSC" & FractionOfTotal>=1)$Gene,
subset(Contributions_InVivo,CellType=="MSC" & FractionOfTotal>=1)$Gene
))


sort(intersect(
  subset(Contributions_InVitro,CellType=="MSC.Diff" & FractionOfTotal>=1)$Gene,
  subset(Contributions_InVivo,CellType=="MSC.Diff" & FractionOfTotal>=1)$Gene
))


sort(intersect(
  subset(Contributions_InVitro,CellType=="VSMC" & FractionOfTotal>=1)$Gene,
  subset(Contributions_InVivo,CellType=="VSMC" & FractionOfTotal>=1)$Gene
))

sort(intersect(
  subset(Contributions_InVitro,CellType=="VEC" & FractionOfTotal>=1)$Gene,
  subset(Contributions_InVivo,CellType=="VEC" & FractionOfTotal>=1)$Gene
))

sort(intersect(
  subset(Contributions_InVitro,CellType=="B cell" & FractionOfTotal>=1)$Gene,
  subset(Contributions_InVivo,CellType=="B cell" & FractionOfTotal>=1)$Gene
))


sort(intersect(
  subset(Contributions_InVitro,CellType=="Macrophage" & FractionOfTotal>=1)$Gene,
  subset(Contributions_InVivo,CellType=="Macrophage" & FractionOfTotal>=1)$Gene
))

sort(intersect(
  subset(Contributions_InVitro,CellType=="T cell" & FractionOfTotal>=1)$Gene,
  subset(Contributions_InVivo,CellType=="T cell" & FractionOfTotal>=1)$Gene
))




### Nnat Network
GOI<-c("Nnat","Agmo","Rtn4","Fbxo45","Adipoq","Atf4","Slc1a3","Thd","Actr3","Arpc2","Dlc1","Pde3b","Met","Car3")

Contributions_InVitro.GOI<-subset(Contributions_InVitro, Gene%in% GOI)

Contributions_InVivo.GOI<-subset(Contributions_InVivo, Gene%in% GOI)



# extract the legend from one of the plots
legend <- get_legend(
  # create some space to the left of the legend
  ggplot( subset(Contributions_InVitro.GOI,CellType %in% c("MSC","MSC.Diff")) ,aes(y=Gene,fill=FractionOfTotal,x=CellType))+
    geom_tile() + theme(legend.box.margin = margin(0, 20, 0, 5))+
    scale_fill_viridis(limits=c(0,1),option="viridis")+
    labs(fill="Fraction of\nexpression")+
    theme(legend.position = "bottom", legend.title = element_text(margin = margin(0,7,0,0)) )
)




InVitro.Adip<-ggplot( subset(Contributions_InVitro.GOI,CellType %in% c("MSC","MSC.Diff")) ,aes(y=Gene,fill=FractionOfTotal,x=CellType))+
  geom_tile()+
  scale_x_discrete(expand = c(0, 0))+
  scale_fill_viridis(limits=c(0,1),option="viridis")+
  scale_y_discrete(expand = c(0, 0))+
  theme(legend.position = "none", plot.margin = unit(c(0.1,0,0.05,0.1), "cm"))+
  xlab("Adipocytes")+
  ggtitle("")


InVitro.Vasc<-ggplot( subset(Contributions_InVitro.GOI,CellType %in% c("VSMC","VEC")) ,aes(y=Gene,fill=FractionOfTotal,x=CellType))+
  geom_tile()+
  scale_fill_viridis(limits=c(0,1),option="viridis")+
  scale_x_discrete(expand = c(0, 0))+
  scale_y_discrete(expand = c(0, 0))+
  theme(legend.position = "none", plot.margin = unit(c(0.1,0,0.05,0), "cm"), axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank())+
  xlab("Vascular")+
  ggtitle("")

InVitro.Immune<-ggplot( subset(Contributions_InVitro.GOI,CellType %in% c("B cell","Macrophage","T cell")) ,aes(y=Gene,fill=FractionOfTotal,x=CellType))+
  geom_tile()+
  scale_fill_viridis(limits=c(0,1),option="viridis")+
  scale_x_discrete(expand = c(0, 0))+
  scale_y_discrete(expand = c(0, 0))+
  theme(legend.position = "none", plot.margin = unit(c(0.1,0,0.05,0), "cm"), axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank())+
  xlab("Immune")+
  ggtitle("")

Plot.InVitro<-plot_grid(Pie.InVitro,InVitro.Adip,InVitro.Vasc,InVitro.Immune, NULL, nrow=1, rel_widths = c(2,2.7,2,3,0.1))






InVivo.Adip<-ggplot( subset(Contributions_InVivo.GOI,CellType %in% c("MSC","MSC.Diff")) ,aes(y=Gene,fill=FractionOfTotal,x=CellType))+
  geom_tile()+
  scale_x_discrete(expand = c(0, 0))+
  scale_fill_viridis(limits=c(0,1),option="viridis")+
  scale_y_discrete(expand = c(0, 0))+
  theme(legend.position = "none", plot.margin = unit(c(0.1,0,0.05,0.1), "cm"))+
  xlab("Adipocytes")+
  ggtitle("")

InVivo.Vasc<-ggplot( subset(Contributions_InVivo.GOI,CellType %in% c("VSMC","VEC")) ,aes(y=Gene,fill=FractionOfTotal,x=CellType))+
  geom_tile()+
  scale_fill_viridis(limits=c(0,1),option="viridis")+
  scale_x_discrete(expand = c(0, 0))+
  scale_y_discrete(expand = c(0, 0))+
  theme(legend.position = "none", plot.margin = unit(c(0.1,0,0.05,0), "cm"), axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank())+
  xlab("Vascular")+
  ggtitle("")

InVivo.Immune<-ggplot( subset(Contributions_InVivo.GOI,CellType %in% c("B cell","Macrophage","T cell")) ,aes(y=Gene,fill=FractionOfTotal,x=CellType))+
  geom_tile()+
  scale_fill_viridis(limits=c(0,1),option="viridis")+
  scale_x_discrete(expand = c(0, 0))+
  scale_y_discrete(expand = c(0, 0))+
  theme(legend.position = "none", plot.margin = unit(c(0.1,0,0.05,0), "cm"), axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank())+
  xlab("Immune")+
  ggtitle("")

Plot.InVivo<-plot_grid(Pie.InVivo,InVivo.Adip,InVivo.Vasc,InVivo.Immune,NULL,nrow=1,rel_widths = c(2,2.7,2,3,0.1))


Plot.Legends<-plot_grid(Pie.Legend,NULL,legend,rel_widths = c(2,0.1,1.25),nrow=1)





Left.Title<-ggplot()+
  annotate("text",x=1,y=1,size=6,label="Cell type composition")+
  theme_void()

Right.Title<-ggplot()+
  annotate("text",x=1,y=1,size=6,label="Cell type contributions")+
  theme_void()

Plot.Titles<-plot_grid(Left.Title,Right.Title,rel_widths = c(1,2),nrow=1)




Full.Panel<-plot_grid(Plot.Titles, Plot.InVitro, Plot.InVivo, Plot.Legends, nrow = 4, rel_heights = c(0.2,1,1,0.45))

ggsave(filename = "CellContributions_Panel_NnatNet.png",plot = Full.Panel,units = "in", height = 5.2, width = 7)







#### Msr1 Net



GOI<-c("Msr1","Cd93","Actr2","P2rx7","Clec10a","Pik3r5","Tpd52","Lrmp","Ly86","Cd84","Arhgap30","Gpnmb","Ccr5","Alox5ap","Rab8b","Myo1f","Cerk","P2ry6","Lgals3","Gas2l3","Naip2")

Contributions_InVitro.GOI<-subset(Contributions_InVitro, Gene%in% GOI)

Contributions_InVivo.GOI<-subset(Contributions_InVivo, Gene%in% GOI)



# extract the legend from one of the plots
legend <- get_legend(
  # create some space to the left of the legend
  ggplot( subset(Contributions_InVitro.GOI,CellType %in% c("MSC","MSC.Diff")) ,aes(y=Gene,fill=FractionOfTotal,x=CellType))+
    geom_tile() + theme(legend.box.margin = margin(0, 20, 0, 5))+
    scale_fill_viridis(limits=c(0,1),option="viridis")+
    labs(fill="Fraction of\nexpression")+
    theme(legend.position = "bottom", legend.title = element_text(margin = margin(0,7,0,0)) )
)




InVitro.Adip<-ggplot( subset(Contributions_InVitro.GOI,CellType %in% c("MSC","MSC.Diff")) ,aes(y=Gene,fill=FractionOfTotal,x=CellType))+
  geom_tile()+
  scale_x_discrete(expand = c(0, 0))+
  scale_fill_viridis(limits=c(0,1),option="viridis")+
  scale_y_discrete(expand = c(0, 0))+
  theme(legend.position = "none", plot.margin = unit(c(0.1,0,0.05,0.1), "cm"))+
  xlab("Adipocytes")+
  ggtitle("")


InVitro.Vasc<-ggplot( subset(Contributions_InVitro.GOI,CellType %in% c("VSMC","VEC")) ,aes(y=Gene,fill=FractionOfTotal,x=CellType))+
  geom_tile()+
  scale_fill_viridis(limits=c(0,1),option="viridis")+
  scale_x_discrete(expand = c(0, 0))+
  scale_y_discrete(expand = c(0, 0))+
  theme(legend.position = "none", plot.margin = unit(c(0.1,0,0.05,0), "cm"), axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank())+
  xlab("Vascular")+
  ggtitle("")

InVitro.Immune<-ggplot( subset(Contributions_InVitro.GOI,CellType %in% c("B cell","Macrophage","T cell")) ,aes(y=Gene,fill=FractionOfTotal,x=CellType))+
  geom_tile()+
  scale_fill_viridis(limits=c(0,1),option="viridis")+
  scale_x_discrete(expand = c(0, 0))+
  scale_y_discrete(expand = c(0, 0))+
  theme(legend.position = "none", plot.margin = unit(c(0.1,0,0.05,0), "cm"), axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank())+
  xlab("Immune")+
  ggtitle("")

Plot.InVitro<-plot_grid(Pie.InVitro,InVitro.Adip,InVitro.Vasc,InVitro.Immune, NULL, nrow=1, rel_widths = c(2,2.7,2,3,0.1))






InVivo.Adip<-ggplot( subset(Contributions_InVivo.GOI,CellType %in% c("MSC","MSC.Diff")) ,aes(y=Gene,fill=FractionOfTotal,x=CellType))+
  geom_tile()+
  scale_x_discrete(expand = c(0, 0))+
  scale_fill_viridis(limits=c(0,1),option="viridis")+
  scale_y_discrete(expand = c(0, 0))+
  theme(legend.position = "none", plot.margin = unit(c(0.1,0,0.05,0.1), "cm"))+
  xlab("Adipocytes")+
  ggtitle("")

InVivo.Vasc<-ggplot( subset(Contributions_InVivo.GOI,CellType %in% c("VSMC","VEC")) ,aes(y=Gene,fill=FractionOfTotal,x=CellType))+
  geom_tile()+
  scale_fill_viridis(limits=c(0,1),option="viridis")+
  scale_x_discrete(expand = c(0, 0))+
  scale_y_discrete(expand = c(0, 0))+
  theme(legend.position = "none", plot.margin = unit(c(0.1,0,0.05,0), "cm"), axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank())+
  xlab("Vascular")+
  ggtitle("")

InVivo.Immune<-ggplot( subset(Contributions_InVivo.GOI,CellType %in% c("B cell","Macrophage","T cell")) ,aes(y=Gene,fill=FractionOfTotal,x=CellType))+
  geom_tile()+
  scale_fill_viridis(limits=c(0,1),option="viridis")+
  scale_x_discrete(expand = c(0, 0))+
  scale_y_discrete(expand = c(0, 0))+
  theme(legend.position = "none", plot.margin = unit(c(0.1,0,0.05,0), "cm"), axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank())+
  xlab("Immune")+
  ggtitle("")

Plot.InVivo<-plot_grid(Pie.InVivo,InVivo.Adip,InVivo.Vasc,InVivo.Immune,NULL,nrow=1,rel_widths = c(2,2.7,2,3,0.1))


Plot.Legends<-plot_grid(Pie.Legend,NULL,legend,rel_widths = c(2,0.1,1.25),nrow=1)





Left.Title<-ggplot()+
  annotate("text",x=1,y=1,size=6,label="Cell type composition")+
  theme_void()

Right.Title<-ggplot()+
  annotate("text",x=1,y=1,size=6,label="Cell type contributions")+
  theme_void()

Plot.Titles<-plot_grid(Left.Title,Right.Title,rel_widths = c(1,2),nrow=1)




Full.Panel<-plot_grid(Plot.Titles, Plot.InVitro, Plot.InVivo, Plot.Legends, nrow = 4, rel_heights = c(0.2,1,1,0.45))

ggsave(filename = "CellContributions_Panel_Msr1Net.png",plot = Full.Panel,units = "in", height = 7.5, width = 7)




#### Plot UMAP
library("Seurat")

data.seurat<-readRDS("data.seurat")




Examples<-c("Sfrp4","Reep6","Myl7","Pde2a","Ly6d","Csf3r","Cd28")


Contributions_InVitro.Examples<-subset(Contributions_InVitro, Gene%in% Examples)

Contributions_InVivo.Examples<-subset(Contributions_InVivo, Gene%in% Examples)



# extract the legend from one of the plots
legend <- get_legend(
  # create some space to the left of the legend
  ggplot( subset(Contributions_InVitro.Examples,CellType %in% c("MSC","MSC.Diff")) ,aes(y=Gene,fill=FractionOfTotal,x=CellType))+
    geom_tile() + theme(legend.box.margin = margin(0, 20, 0, 5))+
    scale_fill_viridis(limits=c(0,1),option="viridis")+
    labs(fill="Fraction of\nexpression")+
    theme(legend.position = "bottom", legend.title = element_text(margin = margin(0,7,0,0)) )
)




InVitro.Adip.Examples<-ggplot( subset(Contributions_InVitro.Examples,CellType %in% c("MSC","MSC.Diff")) ,aes(y=Gene,fill=FractionOfTotal,x=CellType))+
  geom_tile()+
  scale_x_discrete(expand = c(0, 0))+
  scale_fill_viridis(limits=c(0,1),option="viridis")+
  scale_y_discrete(expand = c(0, 0))+
  theme(legend.position = "none", plot.margin = unit(c(0.1,0,0.05,0.1), "cm"))+
  xlab("Adipocytes")+
  ggtitle("")


InVitro.Vasc.Examples<-ggplot( subset(Contributions_InVitro.Examples,CellType %in% c("VSMC","VEC")) ,aes(y=Gene,fill=FractionOfTotal,x=CellType))+
  geom_tile()+
  scale_fill_viridis(limits=c(0,1),option="viridis")+
  scale_x_discrete(expand = c(0, 0))+
  scale_y_discrete(expand = c(0, 0))+
  theme(legend.position = "none", plot.margin = unit(c(0.1,0,0.05,0), "cm"), axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank())+
  xlab("Vascular")+
  ggtitle("")

InVitro.Immune.Examples<-ggplot( subset(Contributions_InVitro.Examples,CellType %in% c("B cell","Macrophage","T cell")) ,aes(y=Gene,fill=FractionOfTotal,x=CellType))+
  geom_tile()+
  scale_fill_viridis(limits=c(0,1),option="viridis")+
  scale_x_discrete(expand = c(0, 0))+
  scale_y_discrete(expand = c(0, 0))+
  theme(legend.position = "none", plot.margin = unit(c(0.1,0,0.05,0), "cm"), axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank())+
  xlab("Immune")+
  ggtitle("")

Plot.InVitro.Examples<-plot_grid(Pie.InVitro,InVitro.Adip.Examples,InVitro.Vasc.Examples,InVitro.Immune.Examples, NULL, nrow=1, rel_widths = c(2,2.7,2,3,0.1))






InVivo.Adip.Examples<-ggplot( subset(Contributions_InVivo.Examples,CellType %in% c("MSC","MSC.Diff")) ,aes(y=Gene,fill=FractionOfTotal,x=CellType))+
  geom_tile()+
  scale_x_discrete(expand = c(0, 0))+
  scale_fill_viridis(limits=c(0,1),option="viridis")+
  scale_y_discrete(expand = c(0, 0))+
  theme(legend.position = "none", plot.margin = unit(c(0.1,0,0.05,0.1), "cm"))+
  xlab("Adipocytes")+
  ggtitle("")

InVivo.Vasc.Examples<-ggplot( subset(Contributions_InVivo.Examples,CellType %in% c("VSMC","VEC")) ,aes(y=Gene,fill=FractionOfTotal,x=CellType))+
  geom_tile()+
  scale_fill_viridis(limits=c(0,1),option="viridis")+
  scale_x_discrete(expand = c(0, 0))+
  scale_y_discrete(expand = c(0, 0))+
  theme(legend.position = "none", plot.margin = unit(c(0.1,0,0.05,0), "cm"), axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank())+
  xlab("Vascular")+
  ggtitle("")

InVivo.Immune.Examples<-ggplot( subset(Contributions_InVivo.Examples,CellType %in% c("B cell","Macrophage","T cell")) ,aes(y=Gene,fill=FractionOfTotal,x=CellType))+
  geom_tile()+
  scale_fill_viridis(limits=c(0,1),option="viridis")+
  scale_x_discrete(expand = c(0, 0))+
  scale_y_discrete(expand = c(0, 0))+
  theme(legend.position = "none", plot.margin = unit(c(0.1,0,0.05,0), "cm"), axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank())+
  xlab("Immune")+
  ggtitle("")

Plot.InVivo.Examples<-plot_grid(Pie.InVivo,InVivo.Adip.Examples,InVivo.Vasc.Examples,InVivo.Immune.Examples,NULL,nrow=1,rel_widths = c(2,2.7,2,3,0.1))


Plot.Legends.Examples<-plot_grid(Pie.Legend,NULL,legend,rel_widths = c(2,0.1,1.25),nrow=1)





Left.Title.Examples<-ggplot()+
  annotate("text",x=1,y=1,size=6,label="Cell type composition")+
  theme_void()

Right.Title.Examples<-ggplot()+
  annotate("text",x=1,y=1,size=6,label="Cell type contributions")+
  theme_void()

Plot.Titles.Examples<-plot_grid(Left.Title.Examples,Right.Title.Examples,rel_widths = c(1,2),nrow=1)




Full.Panel.Examples<-plot_grid(Plot.Titles.Examples, Plot.InVitro.Examples, Plot.InVivo.Examples, Plot.Legends.Examples, nrow = 4, rel_heights = c(0.2,1,1,0.45))

ggsave(filename = "CellContributions_Panel.Examples.png",plot = Full.Panel.Examples,units = "in", height = 4, width = 7)


SpecificGenesExamplesUMAP<-FeaturePlot(data.seurat,Examples, order = TRUE,ncol = 4)+UMAPPlot(data.seurat)

ggsave(filename = "SpecificGenesExamplesUMAP.png",plot = SpecificGenesExamplesUMAP,units = "in", width = 11, height = 4.5)






