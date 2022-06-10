rm(list=ls())
set.seed(0)

require("edgeR")

#Reading counts in

Data<-read.delim("WAT_MergedCounts.tsv",sep=" ",header=TRUE)


Data<-Data[,colnames(Data) != "CCGGACC"]


geneData<-as.matrix(Data[5:nrow(Data),2:ncol(Data)])
rownames(geneData)<-as.character(Data[5:nrow(Data),1])

Info<-read.delim("Sample_Info.txt",sep="\t",header=TRUE)

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
write.table(WATCounts, file='WATCounts.txt', quote=FALSE, sep='\t', row.names=TRUE)


NumerofWATGenes<-dim(WAT_fERdf$counts)[1]

##WAT

WAT_Cross <- topTags(glmLRT(WATfit, coef=2),n=NumerofWATGenes)
WAT_Sex <- topTags(glmLRT(WATfit, coef=4),n=NumerofWATGenes)
WAT_Diet <- topTags(glmLRT(WATfit, coef=3),n=NumerofWATGenes)
WAT_CrossDiet <- topTags(glmLRT(WATfit, coef=c(2,5)),n=NumerofWATGenes)
WAT_CrossSex <- topTags(glmLRT(WATfit, coef=c(2,6)),n=NumerofWATGenes)
WAT_CrossDietSex <- topTags(glmLRT(WATfit, coef=c(2,5,6,8)),n=NumerofWATGenes)

###############




#### Permutation Models
set.seed(0)

iterations<-1:30
PermQuantiles<-NULL
permMeanQuantiles<-NULL
permMeanDiffQuantiles<-NULL

pWAT_Cross <- NULL
pWAT_Sex <- NULL
pWAT_Diet <- NULL
pWAT_CrossDiet <- NULL
pWAT_CrossSex <- NULL
pWAT_CrossDietSex <- NULL



for(i in iterations){
indexes<-sample(1:32,32,replace=FALSE)
pCross<-Info$Cross[BarcodePos][indexes]
pDiet<-Info$Diet[BarcodePos][indexes]
pSex<-Info$Sex[BarcodePos][indexes]

pWAT_fERdf<-DGEList(counts=geneData)
pWATkeep<-rowSums(cpm(pWAT_fERdf)>minCPM) >=minSample ### Maybe change to 4
pWAT_fERdf<-pWAT_fERdf[pWATkeep,keep.lib.sizes=FALSE]
#Normalization
pWAT_fERdf<-calcNormFactors(pWAT_fERdf)
#GLM
pWATdesign <- model.matrix(~pCross*pDiet*pSex)
pWAT_fERdf<-estimateDisp(pWAT_fERdf,pWATdesign)

pWATfit<-glmFit(pWAT_fERdf, pWATdesign)

##pWAT

pWAT_Cross <- rbind(pWAT_Cross,topTags(glmLRT(pWATfit, coef=2),n=NumerofWATGenes)$table)
pWAT_Sex <- rbind(pWAT_Sex,topTags(glmLRT(pWATfit, coef=4),n=NumerofWATGenes)$table)
pWAT_Diet <- rbind(pWAT_Diet,topTags(glmLRT(pWATfit, coef=3),n=NumerofWATGenes)$table)
pWAT_CrossDiet <- rbind(pWAT_CrossDiet,topTags(glmLRT(pWATfit, coef=c(2,5)),n=NumerofWATGenes)$table)
pWAT_CrossSex <- rbind(pWAT_CrossSex,topTags(glmLRT(pWATfit, coef=c(2,6)),n=NumerofWATGenes)$table)
pWAT_CrossDietSex <- rbind(pWAT_CrossDietSex,topTags(glmLRT(pWATfit, coef=c(2,5,6,8)),n=NumerofWATGenes)$table)

###############

PermQuantiles<-rbind(PermQuantiles,quantile(pWAT_CrossDietSex$LR,type=4,probs=seq(0,1,0.01)))
permMeanQuantiles<-rbind(permMeanQuantiles,apply(PermQuantiles,2,mean))
permMeanDiffQuantiles<-rbind(permMeanDiffQuantiles,permMeanQuantiles[i-1,]-permMeanQuantiles[i,])

}

tiff(filename="WAT_Permutation_mean_deviation.tiff",width = 600, height = 600, res=100)
plot(y=permMeanDiffQuantiles[,51],iterations[1:length(iterations)-1],type='l',ylab="Quantile Deviation from Mean",xlab="Iteration")
dev.off()


library(lsr)
NullDist_Cross<-ecdf(pWAT_Cross$LR)
NullDist_Sex<-ecdf(pWAT_Sex$LR)
NullDist_Diet<-ecdf(pWAT_Diet$LR)

NullDist_CrossDiet<-ecdf(pWAT_CrossDiet$LR)
NullDist_CrossSex<-ecdf(pWAT_CrossSex$LR)

NullDist_CrossDietSex<-ecdf(pWAT_CrossDietSex$LR)


WAT_Cross_P<-unlist(lapply(abs(WAT_Cross$table$LR),function(t) 1-NullDist_Cross(t)))
WAT_Sex_P<-unlist(lapply(abs(WAT_Sex$table$LR),function(t) 1-NullDist_Cross(t)))
WAT_Diet_P<-unlist(lapply(abs(WAT_Diet$table$LR),function(t) 1-NullDist_Cross(t)))
WAT_CrossDiet_P<-unlist(lapply(abs(WAT_CrossDiet$table$LR),function(t) 1-NullDist_Cross(t)))
WAT_CrossSex_P<-unlist(lapply(abs(WAT_CrossSex$table$LR),function(t) 1-NullDist_Cross(t)))
WAT_CrossDietSex_P<-unlist(lapply(abs(WAT_CrossDietSex$table$LR),function(t) 1-NullDist_Cross(t)))


png(filename="LR_Distributions.png", width=500,height=500,type=c("cairo"))

par(mfrow=c(2,2))

plot(density(as.numeric(as.character(pWAT_Cross$LR
)), na.rm=TRUE),col="orange",lwd=6,main="LR-statistic Distributions")

lines(density(as.numeric(as.character(WAT_Cross$table$LR)), na.rm=TRUE),col="blue",lwd=2)


plot(density(as.numeric(as.character(pWAT_CrossSex$LR
)), na.rm=TRUE),col="orange",lwd=6,main="LR-statistic Distributions")

lines(density(as.numeric(as.character(WAT_CrossSex$table$LR)), na.rm=TRUE),col="blue",lwd=2)


plot(density(as.numeric(as.character(pWAT_CrossDiet$LR
)), na.rm=TRUE),col="orange",lwd=6,main="LR-statistic Distributions")

lines(density(as.numeric(as.character(WAT_CrossDiet$table$LR)), na.rm=TRUE),col="blue",lwd=2)


plot(density(as.numeric(as.character(pWAT_CrossDietSex$LR)), na.rm=TRUE),col="orange",lwd=6,main="LR-statistic Distributions")

lines(density(as.numeric(as.character(WAT_CrossDietSex$table$LR)), na.rm=TRUE),col="blue",lwd=2)


dev.off()



WAT_Cross$table$P<-WAT_Cross_P
WAT_Sex$table$P<-WAT_Sex_P
WAT_Diet$table$P<-WAT_Diet_P
WAT_CrossDiet$table$P<-WAT_CrossDiet_P
WAT_CrossSex$table$P<-WAT_CrossSex_P
WAT_CrossDietSex$table$P<-WAT_CrossDietSex_P

minFC<-1

passed_WAT_Cross <- WAT_Cross$table[WAT_Cross_P<=0.05 & abs(WAT_Cross$table$logFC)>=minFC,]
passed_WAT_Sex <- WAT_Sex$table[WAT_Sex_P<=0.05 & abs(WAT_Sex$table$logFC)>=minFC,]
passed_WAT_Diet <- WAT_Diet$table[WAT_Diet_P<=0.05 & abs(WAT_Diet$table$logFC)>=minFC,]

passed_WAT_CrossDiet <- WAT_CrossDiet$table[WAT_CrossDiet_P<=0.05 & (abs(WAT_CrossDiet$table$logFC.CrossSXL)>=minFC | abs(WAT_CrossDiet$table$logFC.CrossSXL.DietL)>=minFC ),]

passed_WAT_CrossSex<- WAT_CrossSex$table[WAT_CrossSex_P<=0.05 & (abs(WAT_CrossSex$table$logFC.CrossSXL)>=minFC | abs(WAT_CrossSex$table$logFC.CrossSXL.SexM)>=minFC),]

passed_WAT_CrossDietSex<- WAT_CrossDietSex$table[WAT_CrossDietSex_P<=0.05 & (abs(WAT_CrossDietSex$table$logFC.CrossSXL)>=minFC | abs(WAT_CrossDietSex$table$logFC.CrossSXL.SexM)>=minFC | abs(WAT_CrossDietSex $table$logFC.CrossSXL.DietL)>=minFC | abs(WAT_CrossDietSex$table$logFC.CrossSXL.DietL.SexM)>=minFC),]


length(unique(c(
rownames(passed_WAT_Cross),
rownames(passed_WAT_CrossDiet),
rownames(passed_WAT_CrossSex),
rownames(passed_WAT_CrossDietSex)
)))


NumerofWATGenes


write.table(passed_WAT_Cross, file='passed_WAT_Cross.txt', quote=FALSE, sep='\t', row.names=TRUE)
write.table(passed_WAT_Sex, file='passed_WAT_Sex.txt', quote=FALSE, sep='\t', row.names=TRUE)
write.table(passed_WAT_Diet, file='passed_WAT_Diet.txt', quote=FALSE, sep='\t', row.names=TRUE)

write.table(passed_WAT_CrossDiet, file='passed_WAT_CrossDiet.txt', quote=FALSE, sep='\t', row.names=TRUE)
write.table(passed_WAT_CrossSex, file='passed_WAT_CrossSex.txt', quote=FALSE, sep='\t', row.names=TRUE)
write.table(passed_WAT_CrossDietSex, file='passed_WAT_CrossDietSex.txt', quote=FALSE, sep='\t', row.names=TRUE)



write.table(WAT_Cross$table, file='All_WAT_Cross.txt', quote=FALSE, sep='\t', row.names=TRUE)
write.table(WAT_Sex$table, file='All_WAT_Sex.txt', quote=FALSE, sep='\t', row.names=TRUE)
write.table(WAT_Diet$table, file='All_WAT_Diet.txt', quote=FALSE, sep='\t', row.names=TRUE)

write.table(WAT_CrossDiet$table, file='All_WAT_CrossDiet.txt', quote=FALSE, sep='\t', row.names=TRUE)
write.table(WAT_CrossSex$table, file='All_WAT_CrossSex.txt', quote=FALSE, sep='\t', row.names=TRUE)
write.table(WAT_CrossDietSex$table, file='All_WAT_CrossDietSex.txt', quote=FALSE, sep='\t', row.names=TRUE)




write.table(
unique(c(
rownames(passed_WAT_Cross),
rownames(passed_WAT_CrossDiet),
rownames(passed_WAT_CrossSex),
rownames(passed_WAT_CrossDietSex)
))
,file='passed_All.txt', quote=FALSE,col.names=FALSE,row.names=FALSE
)





PostionData<-read.delim("LG_GeneBodyAnnotations.txt",header=FALSE)
colnames(PostionData)<-c("Chr","Start","Stop","Gene")


IDCon<-read.delim("IDconversions.txt")
matchedIDs<-subset(IDCon,toupper(Gene.stable.ID) %in% toupper(rownames(geneData)) & !is.na(NCBI.gene.ID))

#posWATCounts<-WATCounts[rownames(WATCounts) %in% matchedIDs$Gene.stable.ID[matchedIDs$Gene.name %in% PostionData$Gene],]



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

require("ggplot2")
require("cowplot")
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
wat_data_hf<-subset(wat_data,Sex=="F" & Diet=="H")


GENE<-"ENSMUSG00000012187" #Mogat1
GOI_wat_data_hf<-wat_data_hf[,c(1:3,(1:ncol(wat_data_hf))[colnames(wat_data_hf) == GENE])]
colnames(GOI_wat_data_hf)[colnames(GOI_wat_data_hf)==GENE]<-"Expression"
GOIDF<-data_summary(GOI_wat_data_hf, varname="Expression", groupnames=c("Cross","Sex"))
 
tiff(filename="Mogat1.tiff",width = 600, height = 600, res=100)
ggplot(data= GOIDF, aes(x= Cross, y= Expression) ) +
  geom_bar(stat="identity", width=0.5) + 
  geom_errorbar(aes(ymin= Expression-sem, ymax= Expression +sem), width=.2, position=position_dodge(0))
dev.off()


GENE<-"ENSMUSG00000067786" #Nnat
GOI_wat_data_hf<-wat_data_hf[,c(1:3,(1:ncol(wat_data_hf))[colnames(wat_data_hf) == GENE])]
colnames(GOI_wat_data_hf)[colnames(GOI_wat_data_hf)==GENE]<-"Expression"
GOIDF<-data_summary(GOI_wat_data_hf, varname="Expression", groupnames=c("Cross","Sex"))

tiff(filename="Nnat.tiff",width = 600, height = 600, res=100)
ggplot(data= GOIDF, aes(x= Cross, y= Expression)) +
  geom_bar(stat="identity", width=0.5) + 
  geom_errorbar(aes(ymin= Expression-sem, ymax= Expression +sem), width=.2, position=position_dodge(0))
dev.off()



