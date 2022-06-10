QTL_data<-read.table("genes_in_POE_QTL.BED",header=FALSE)
colnames(QTL_data)<-c("Chr","Start","Stop","Gene","QTL")
poeQTLgenes<-unique(QTL_data$Gene)


poefiltEpiData<-subset(filtEpi_FDR, gA %in% poeQTLgenes)

dim(poefiltEpiData)

length(unique(poefiltEpiData$gA))
length(unique(poefiltEpiData$gB))

write.table(poefiltEpiData, file='poefiltEpiData.tsv', quote=FALSE, sep='\t', row.names=FALSE)



