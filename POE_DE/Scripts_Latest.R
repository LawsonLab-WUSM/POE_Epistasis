


rm(list=ls())
set.seed(0)

require("edgeR")
require("cowplot")
require("ggplot2")


COUNT_FILE_NAME<-"WAT_MergedCounts.tsv"
SAMPLE_INFO_FILE_NAME<-"Sample_Info.txt"
NORMALIZED_COUNT_FILE_NAME<-'WATCounts.txt'
NAME_CONVERSION_FILE<-"IDconversions_LGSM.txt"
GENES_IN_QTL_FILE<-"genes_in_POE_QTL.BED"
ANNOTATION_FILE<-"LG_GeneBodyAnnotations.txt"
NCBI_ID_CONVERSION_FILE<-"IDconversions.txt"

LIST_OF_BARCODES_TO_OMIT<-c("CCGGACC")


#POE QTL filter
NameConversion<-read.delim(NAME_CONVERSION_FILE)

QTL_data<-read.table(GENES_IN_QTL_FILE,header=FALSE)


colnames(QTL_data)<-c("Chr","Start","Stop","Gene","QTL")
QTL_data$ENSEMBL<-NameConversion$Gene.stable.ID[match(QTL_data$Gene, NameConversion$Gene.name)]
poeQTLgenes<-as.character(unique(QTL_data$ENSEMBL))


PostionData<-read.delim(ANNOTATION_FILE,header=FALSE)
colnames(PostionData)<-c("Chr","Start","Stop","Gene")


IDCon<-read.delim(NCBI_ID_CONVERSION_FILE)



#Reading counts in

Data<-read.delim(COUNT_FILE_NAME,sep=" ",header=TRUE)


Data<-Data[,!(colnames(Data) %in% LIST_OF_BARCODES_TO_OMIT)]


geneData<-as.matrix(Data[5:nrow(Data),2:ncol(Data)])
rownames(geneData)<-as.character(Data[5:nrow(Data),1])


Info<-read.delim(SAMPLE_INFO_FILE_NAME,sep="\t",header=TRUE)

BarcodePos<-match(colnames(geneData),Info$GTAC.Tag)

Cross<-Info$Cross[BarcodePos]
Diet<-Info$Diet[BarcodePos]
Sex<-Info$Sex[BarcodePos]

WAT_fERdf<-DGEList(counts=geneData)

#Filtering

minCPM<-10
minSample<-16

WATkeep<-rowSums(cpm(WAT_fERdf)>minCPM) >=minSample ### Maybe change to 4


WAT_fERdf<-WAT_fERdf[WATkeep,keep.lib.sizes=FALSE]

#Normalization

WAT_fERdf<-calcNormFactors(WAT_fERdf)

#GLM

WATdesign <- model.matrix(~Cross*Diet*Sex)

WAT_fERdf<-estimateDisp(WAT_fERdf,WATdesign)

WATfit<-glmFit(WAT_fERdf, WATdesign)
WATCounts<-cpm(WAT_fERdf,normalized=T)
write.table(WATCounts, file=NORMALIZED_COUNT_FILE_NAME, quote=FALSE, sep='\t', row.names=TRUE)





NumerofWATGenes<-dim(WAT_fERdf$counts)[1]

##WAT

WAT_Cross <- topTags(glmLRT(WATfit, coef=2),n=NumerofWATGenes)
WAT_Sex <- topTags(glmLRT(WATfit, coef=4),n=NumerofWATGenes)
WAT_Diet <- topTags(glmLRT(WATfit, coef=3),n=NumerofWATGenes)
WAT_CrossDiet <- topTags(glmLRT(WATfit, coef=c(5)),n=NumerofWATGenes)
WAT_CrossSex <- topTags(glmLRT(WATfit, coef=c(6)),n=NumerofWATGenes)
WAT_CrossDietSex <- topTags(glmLRT(WATfit, coef=c(8)),n=NumerofWATGenes)
WAT_CrossFull <- topTags(glmLRT(WATfit, coef=c(2,5,6,8)),n=NumerofWATGenes)




###############

SignificanceThreshold <-0.1

subset(WAT_Cross$table, FDR<= SignificanceThreshold)

subset(WAT_CrossDiet$table, FDR<= SignificanceThreshold)

subset(WAT_CrossSex$table, FDR<= SignificanceThreshold)

subset(WAT_CrossDietSex$table, FDR<= SignificanceThreshold)

subset(WAT_CrossFull$table, FDR<= SignificanceThreshold)


library(qvalue)



nrow(subset(WAT_Cross$table, FDR<= SignificanceThreshold))

nrow(subset(WAT_CrossDiet$table, FDR<= SignificanceThreshold))

nrow(subset(WAT_CrossSex$table, FDR<= SignificanceThreshold))

nrow(subset(WAT_CrossDietSex$table, FDR<= SignificanceThreshold))

nrow(subset(WAT_CrossFull$table, FDR<= SignificanceThreshold))



Qval_Sex<-qvalue(WAT_Sex$table$PValue)
Qval_Diet<-qvalue(WAT_Diet$table$PValue)

Qval_Cross<-qvalue(WAT_Cross$table$PValue)
Qval_CrossDiet<-qvalue(WAT_CrossDiet$table$PValue)
Qval_CrossSex<-qvalue(WAT_CrossSex$table$PValue)
Qval_CrossDietSex<-qvalue(WAT_CrossDietSex$table$PValue)
Qval_CrossFull<-qvalue(WAT_CrossFull$table$PValue)


png(filename="Pval_Dist_Sex.png", width=500,height=500,type=c("cairo"))
hist(Qval_Sex)
dev.off()

png(filename="Pval_Dist_Diet.png", width=500,height=500,type=c("cairo"))
hist(Qval_Diet)
dev.off()


png(filename="Pval_Dist_Cross.png", width=500,height=500,type=c("cairo"))
hist(Qval_Cross)
dev.off()

png(filename="Pval_Dist_CrossDiet.png", width=500,height=500,type=c("cairo"))
hist(Qval_CrossDiet)
dev.off()


png(filename="Pval_Dist_CrossSex.png", width=500,height=500,type=c("cairo"))
hist(Qval_CrossSex)
dev.off()

png(filename="Pval_Dist_CrossDietSex.png", width=500,height=500,type=c("cairo"))
hist(Qval_CrossDietSex)
dev.off()


png(filename="Pval_Dist_CrossFull.png", width=500,height=500,type=c("cairo"))
hist(Qval_CrossFull)
dev.off()



png(filename="MultipleTestsCorrection_Sex.png", width=500,height=500,type=c("cairo"))
plot(Qval_Sex)
dev.off()

png(filename="MultipleTestsCorrection_Diet.png", width=500,height=500,type=c("cairo"))
plot(Qval_Diet)
dev.off()

png(filename="MultipleTestsCorrection_Cross.png", width=500,height=500,type=c("cairo"))
plot(Qval_Cross)
dev.off()


png(filename="MultipleTestsCorrection_CrossDiet.png", width=500,height=500,type=c("cairo"))
plot(Qval_CrossDiet)
dev.off()


png(filename="MultipleTestsCorrection_CrossSex.png", width=500,height=500,type=c("cairo"))
plot(Qval_CrossSex)
dev.off()


png(filename="MultipleTestsCorrection_CrossDietSex.png", width=500,height=500,type=c("cairo"))
plot(Qval_CrossDietSex)
dev.off()

png(filename="MultipleTestsCorrection_CrossFull.png", width=500,height=500,type=c("cairo"))
plot(Qval_CrossFull)
dev.off()




passedSex<-rownames(WAT_Sex)[Qval_Sex$qvalues<0.1]
passedDiet<-rownames(WAT_Diet)[Qval_Diet$qvalues<0.1]


passedCross<-rownames(WAT_Cross)[Qval_Cross$qvalues<0.1]

passedCrossDiet<-rownames(WAT_CrossDiet)[qvalue(WAT_CrossDiet $table$PValue)$qvalues<0.1]

passedCrossSex<-rownames(WAT_CrossSex)[qvalue(WAT_CrossSex $table$PValue)$qvalues<0.1]

passedCrossDietSex<-rownames(WAT_CrossDietSex)[qvalue(WAT_CrossDietSex $table$PValue)$qvalues<0.1]

passedFull<-rownames(WAT_CrossFull)[qvalue(WAT_CrossFull $table$PValue)$qvalues<0.1]


passedAll<-unique(c(passedCross, passedCrossDiet, passedCrossSex, passedCrossDietSex, passedFull))


passedAll[passedAll %in% poeQTLgenes]


length(passedAll)
sum(passedAll %in% passedSex)
sum(passedAll %in% passedDiet)


WAT_Cross$table$Q<-Qval_Cross$qvalues
WAT_Sex$table$Q<-Qval_Sex$qvalues
WAT_Diet$table$Q<-Qval_Diet$qvalues
WAT_CrossDiet$table$Q<-Qval_CrossDiet$qvalues
WAT_CrossSex$table$Q<-Qval_CrossSex$qvalues
WAT_CrossDietSex$table$Q<-Qval_CrossDietSex$qvalues
WAT_CrossFull$table$Q<-Qval_CrossFull$qvalues

minFC<-log2(1.5)




Volcano_Diet<-ggplot(WAT_Diet$table,aes(y=-log10(Q),x=logFC))+
geom_point(size=0.5)+
geom_hline(yintercept= -log10(SignificanceThreshold),col='grey')+
geom_vline(xintercept= minFC, col='grey')+
geom_vline(xintercept= -minFC, col='grey')+
ylab("-Log10( p-Value )")+
xlab("logFC")+
theme(legend.position = "none")+
scale_color_manual(values=c("#999999", "#E69F00"))+
ggtitle("Diet")

ggsave("Volcano_Diet.tiff", Volcano_Diet)





Volcano_Sex<-ggplot(WAT_Sex$table,aes(y=-log10(Q),x=logFC))+
geom_point(size=0.5)+
geom_hline(yintercept= -log10(SignificanceThreshold),col='grey')+
geom_vline(xintercept= minFC, col='grey')+
geom_vline(xintercept= -minFC, col='grey')+
ylab("-Log10( p-Value )")+
xlab("logFC")+
theme(legend.position = "none")+
scale_color_manual(values=c("#999999", "#E69F00"))+
ggtitle("Sex")

ggsave("Volcano_Sex.tiff", Volcano_Sex)






Volcano_Cross<-ggplot(WAT_Cross$table,aes(y=-log10(Q),x=logFC))+
geom_point(size=0.5)+
geom_hline(yintercept= -log10(SignificanceThreshold),col='grey')+
geom_vline(xintercept= minFC, col='grey')+
geom_vline(xintercept= -minFC, col='grey')+
ylab("-Log10( p-Value )")+
xlab("logFC")+
theme(legend.position = "none")+
scale_color_manual(values=c("#999999", "#E69F00"))+
ggtitle("Cross")

ggsave("Volcano_Cross.tiff", Volcano_Cross)



Volcano_CrossDiet<-ggplot(WAT_CrossDiet$table,aes(y=-log10(Q),x=logFC))+
geom_point(size=0.5)+
geom_hline(yintercept= -log10(SignificanceThreshold),col='grey')+
geom_vline(xintercept= minFC, col='grey')+
geom_vline(xintercept= -minFC, col='grey')+
ylab("-Log10( p-Value )")+
xlab("logFC")+
theme(legend.position = "none")+
scale_color_manual(values=c("#999999", "#E69F00"))+
ggtitle("Cross:Diet")

ggsave("Volcano_CrossDiet.tiff", Volcano_CrossDiet)



Volcano_CrossSex<-ggplot(WAT_CrossSex$table,aes(y=-log10(Q),x=logFC))+
geom_point(size=0.5)+
geom_hline(yintercept= -log10(SignificanceThreshold),col='grey')+
geom_vline(xintercept= minFC, col='grey')+
geom_vline(xintercept= -minFC, col='grey')+
ylab("-Log10( p-Value )")+
xlab("logFC")+
theme(legend.position = "none")+
scale_color_manual(values=c("#999999", "#E69F00"))+
ggtitle("Cross:Sex")

ggsave("Volcano_CrossSex.tiff", Volcano_CrossSex)



Volcano_CrossDietSex<-ggplot(WAT_CrossDietSex$table,aes(y=-log10(Q),x=logFC))+
geom_point(size=0.5)+
geom_hline(yintercept= -log10(SignificanceThreshold),col='grey')+
geom_vline(xintercept= minFC, col='grey')+
geom_vline(xintercept= -minFC, col='grey')+
ylab("-Log10( p-Value )")+
xlab("logFC")+
theme(legend.position = "none")+
scale_color_manual(values=c("#999999", "#E69F00"))+
ggtitle("Cross:Diet:Sex")

ggsave("Volcano_CrossDietSex.tiff", Volcano_CrossDietSex)


passed_WAT_Cross <- WAT_Cross$table[WAT_Cross$table$Q<= SignificanceThreshold & abs(WAT_Cross$table$logFC)>=minFC,]
passed_WAT_Sex <- WAT_Sex$table[WAT_Sex$table$Q<= SignificanceThreshold & abs(WAT_Sex$table$logFC)>=minFC,]
passed_WAT_Diet <- WAT_Diet$table[WAT_Diet$table$Q<= SignificanceThreshold & abs(WAT_Diet$table$logFC)>=minFC,]

passed_WAT_CrossDiet <- WAT_CrossDiet$table[WAT_CrossDiet$table$Q<= SignificanceThreshold & (abs(WAT_CrossDiet$table$logFC)>=minFC | abs(WAT_CrossDiet$table$logFC)>=minFC ),]

passed_WAT_CrossSex<- WAT_CrossSex$table[WAT_CrossSex$table$Q<= SignificanceThreshold & (abs(WAT_CrossSex$table$logFC)>=minFC),]

passed_WAT_CrossDietSex<- WAT_CrossDietSex$table[WAT_CrossDietSex$table$Q<= SignificanceThreshold & (abs(WAT_CrossDietSex$table$logFC)>=minFC),]


length(unique(c(
rownames(passed_WAT_Cross),
rownames(passed_WAT_CrossDiet),
rownames(passed_WAT_CrossSex),
rownames(passed_WAT_CrossDietSex)
)))


NumerofWATGenes





write.table(passed_WAT_Cross, file='passed_Cross.txt', quote=FALSE, sep='\t', row.names=TRUE)
write.table(passed_WAT_Sex, file='passed_Sex.txt', quote=FALSE, sep='\t', row.names=TRUE)
write.table(passed_WAT_Diet, file='passed_Diet.txt', quote=FALSE, sep='\t', row.names=TRUE)

write.table(passed_WAT_CrossDiet, file='passed_CrossDiet.txt', quote=FALSE, sep='\t', row.names=TRUE)
write.table(passed_WAT_CrossSex, file='passed_CrossSex.txt', quote=FALSE, sep='\t', row.names=TRUE)
write.table(passed_WAT_CrossDietSex, file='passed_CrossDietSex.txt', quote=FALSE, sep='\t', row.names=TRUE)



write.table(WAT_Cross$table, file='All_Cross.txt', quote=FALSE, sep='\t', row.names=TRUE)
write.table(WAT_Sex$table, file='All_Sex.txt', quote=FALSE, sep='\t', row.names=TRUE)
write.table(WAT_Diet$table, file='All_Diet.txt', quote=FALSE, sep='\t', row.names=TRUE)

write.table(WAT_CrossDiet$table, file='All_CrossDiet.txt', quote=FALSE, sep='\t', row.names=TRUE)
write.table(WAT_CrossSex$table, file='All_CrossSex.txt', quote=FALSE, sep='\t', row.names=TRUE)
write.table(WAT_CrossDietSex$table, file='All_CrossDietSex.txt', quote=FALSE, sep='\t', row.names=TRUE)




passed_All<-unique(c(
rownames(passed_WAT_Cross),
rownames(passed_WAT_CrossDiet),
rownames(passed_WAT_CrossSex),
rownames(passed_WAT_CrossDietSex)
))

write.table(
passed_All,file='passed_All.txt', quote=FALSE,col.names=FALSE,row.names=FALSE
)






GeneName_Passed_All<-as.character(NameConversion$Gene.name[match(passed_All,NameConversion$Gene.stable.ID)])

write.table(
GeneName_Passed_All,file='GeneName_Passed_All.txt', quote=FALSE,col.names=FALSE,row.names=FALSE
)



write.table(
rownames(WAT_CrossDiet),
file='Expressed_All.txt', quote=FALSE,col.names=FALSE,row.names=FALSE
)




matchedIDs<-subset(IDCon,toupper(Gene.stable.ID) %in% toupper(rownames(geneData)) & !is.na(NCBI.gene.ID))






POEscores<-cbind(
rep(0,nrow(WATCounts))
,rep(0,nrow(WATCounts))
,rep(0,nrow(WATCounts))
,rep(0,nrow(WATCounts))
)
colnames(POEscores)<-c("HF","LF","HM","LM")
rownames(POEscores)<-rownames(WATCounts)

Context<-cbind(c("H","L","H","L"),c("F","F","M","M"))


for(i in 1:4){
	POEscores[,i]<-log2(apply(WATCounts[,Diet==Context[i,1] & Sex==Context[i,2] & Cross=="SXL"],1,mean)/apply(WATCounts[,Diet==Context[i,1] & Sex==Context[i,2] & Cross=="LXS"],1,mean))
}

write.table(POEscores, file='POEscores.txt', quote=FALSE, sep='\t', row.names=TRUE)

data_summary <- function(data, varname, groupnames){
  require(plyr)
  summary_func <- function(x, col){
    c(mean = mean(x[[col]], na.rm=TRUE),
      sd = sd(x[[col]], na.rm=TRUE),
      sem = sd(x[[col]], na.rm=TRUE)/(length(x[[col]])^0.5)
      )
  }
  data_sum<-ddply(data, groupnames, .fun=summary_func,
                  varname)
  data_sum <- rename(data_sum, c("mean" = varname))
 return(data_sum)
}

wat_data<-data.frame(Cross,Sex,Diet,t(WATCounts))








DE_Bar_Dir<-"DE_Barplots"
dir.create(DE_Bar_Dir)


for(GENE in passed_All){
	
GOI_wat_data<-wat_data[,c(1:3,(1:ncol(wat_data))[colnames(wat_data) == GENE])]
colnames(GOI_wat_data)[colnames(GOI_wat_data)==GENE]<-"Expression"
GOIDF<-data_summary(GOI_wat_data, varname="Expression", groupnames=c("Cross","Sex","Diet"))


PLOT<-ggplot(data= GOIDF, aes(x= Cross, y= Expression) ) +
  geom_bar(stat="identity", width=0.5) + 
  geom_errorbar(aes(ymin=Expression - sem, ymax=Expression + sem), width=.2, position=position_dodge(0))+
  facet_grid(~Diet*Sex)


ggsave( paste(DE_Bar_Dir,"/DE_Plots_",GENE,".tiff",sep=""), PLOT, units="in", width=6, height=3)


}


