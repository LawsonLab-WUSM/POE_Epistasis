
runTests<-function(){

indexes<-sample(1:32,32,replace=FALSE)
pCross<-Info$Cross[BarcodePos][indexes]
pDiet<-Info$Diet[BarcodePos][indexes]
pSex<-Info$Sex[BarcodePos][indexes]

# pWAT_fERdf<-DGEList(counts=geneData)
pWAT_fERdf<-DGEList(counts=geneData)

pGeneData<-geneData[sample(1:nrow(geneData),nrow(geneData),replace=FALSE),]

rownames(pGeneData)<-rownames(geneData)


pWATkeep<-rowSums(cpm(pWAT_fERdf)>minCPM) >=minSample ### Maybe change to 4
pWAT_fERdf<-pWAT_fERdf[pWATkeep,keep.lib.sizes=FALSE]
#Normalization
pWAT_fERdf<-calcNormFactors(pWAT_fERdf)
#GLM
pWATdesign <- model.matrix(~pCross*pDiet*pSex)
pWAT_fERdf<-estimateDisp(pWAT_fERdf,pWATdesign)

pWATfit<-glmFit(pWAT_fERdf, pWATdesign)


pWAT_Cross <- topTags(glmLRT(pWATfit, coef=2),n=NumerofWATGenes)$table
pWAT_Sex <- topTags(glmLRT(pWATfit, coef=4),n=NumerofWATGenes)$table
pWAT_Diet <- topTags(glmLRT(pWATfit, coef=3),n=NumerofWATGenes)$table
pWAT_CrossDiet <- topTags(glmLRT(pWATfit, coef=c(5)),n=NumerofWATGenes)$table
pWAT_CrossSex <- topTags(glmLRT(pWATfit, coef=c(6)),n=NumerofWATGenes)$table
pWAT_CrossDietSex <- topTags(glmLRT(pWATfit, coef=c(8)),n=NumerofWATGenes)$table
pWAT_CrossFull <- topTags(glmLRT(pWATfit, coef=c(2,5,6,8)),n=NumerofWATGenes)$table



data.frame(
Cross=pWAT_Cross$PValue,
Sex=pWAT_Sex$PValue,
Diet=pWAT_Diet$PValue,
CrossDiet=pWAT_CrossDiet$PValue,
CrossSex=pWAT_CrossSex$PValue,
CrossDietSex=pWAT_CrossDietSex$PValue,
CrossFull=pWAT_CrossFull$PValue
)


}





