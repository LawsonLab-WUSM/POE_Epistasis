
set.seed(0)
rm(list=ls())

library(lsr)
library(qvalue)
library(ggplot2)
library(cowplot)
library(plyr)


#### Put File Names Below
Netmm9_BEDFile<-"WATNet_mm9.BED" ### Transcript Positions Bed file
NetAdditiveScoresFile<-"WATNetAdditiveScores.csv" ### F16 Additive Scores
NetImprintingScoresFile<-"WATNetImprintingScores.csv"	### F16 Imprinting Scores
PhenotypesFile<-"F16Phenotypes.csv"	### F16 Phenotypes
NetTableFile<-"WATNetTable.tsv" ### Relavent DE~ASE gene pairs AND relavent Phenotypes from QTL DE gene falls in
SNPBEDFile<-"F16_SNPS_within_3Mb_of_WATNet.BED" ### List of markers for each transcript

GenesInPOEQTL<-as.character(unique(read.table("genes_in_POE_QTL.BED")[,4]))

#### SNPBEDFile format: [Marker Chr] [Marker Start] [Marker Stop] [Transcript ID] [Gene Name]

#### NetTableFile format: [DE Gene] [ASE Gene] [pvalue] [F1 Trait Name] [F16 Trait Name]

##Import Data
{

Netmm9<-read.delim(Netmm9_BEDFile,header=FALSE,sep="\t")
colnames(Netmm9)<-c("Chr","Start","Stop","ENSEMBL","Gene")

### This is identifying the center point of each gene
NetMedian<-data.frame(matrix(unlist( lapply(1:length(unique(Netmm9$Gene)), function(n) c( mean(c(min(unlist(subset(Netmm9,Gene==unique(Netmm9$Gene)[n])[2:3])),max(unlist(subset(Netmm9,Gene==unique(Netmm9$Gene)[n])[2:3])))), as.character(unique(Netmm9$Gene)[n]) ) )), byrow=TRUE, ncol=2))

colnames(NetMedian)<-c("Locus","Gene")

NetAdditiveScores<-read.csv(NetAdditiveScoresFile)
NetAdditiveScores$LocusID<-apply(NetAdditiveScores[,c(2,5)],1,paste,collapse="_")
NetAdditiveScores$Marker.ID<-as.character(NetAdditiveScores$Marker.ID)

NetImprintingScores<-read.csv(NetImprintingScoresFile)
NetImprintingScores$LocusID<-apply(NetImprintingScores[,c(2,5)],1,paste,collapse="_")
NetImprintingScores$Marker.ID<-as.character(NetImprintingScores$Marker.ID)

F16Phenotypes<-read.csv(PhenotypesFile)

Net<-read.delim(NetTableFile,header=TRUE,sep="\t")

unique(Net[,4:5])


#Net <-subset(Net, Target.Gene %in% GenesInPOEQTL)

#nrow(Net)
#Net<-subset(Net,Trait %in% c("GlucWt1","GlucWt2","CHOL","TG","NecWeight","Insulin","AUC2","AUC1",'T0_1',"T0_2"))



GeneSNPs<-read.delim(SNPBEDFile,header=FALSE,sep="\t")
colnames(GeneSNPs)<-c("Chr","Start","Stop","ENSEMBL","Gene")

GeneSNPs$LocusID<-apply(GeneSNPs[1:2],1,paste,collapse="_")

###### SUPER IMPORTANT!!!!! Make sure R is treating phenotypes as numbers NOT factors
F16Phenotypes$Animal<-unlist(lapply(F16Phenotypes$Animal,function(Z) paste(c("X",Z),collapse="")))

}



### Transforming phenotype
{


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

}


### Covariate Screen
{

head(tF16Phenotypes)
sF16Phenotypes<-tF16Phenotypes

sF16Phenotypes$GlucWt1[as.numeric(rownames(data.frame(residuals(aov(GlucWt1~Sex*Diet*+Litter_Size+Family,sF16Phenotypes)))))]<-residuals(aov(GlucWt1~Sex*Diet*+Litter_Size+Family,sF16Phenotypes))

sF16Phenotypes$GlucWt2[as.numeric(rownames(data.frame(residuals(aov(GlucWt2~Sex*Diet*+Litter_Size+Family,sF16Phenotypes)))))]<-residuals(aov(GlucWt2~Sex*Diet*+Litter_Size+Family,sF16Phenotypes))

sF16Phenotypes$T0_1[as.numeric(rownames(data.frame(residuals(aov(T0_1~Sex*Diet*+Litter_Size+Family,sF16Phenotypes)))))]<-residuals(aov(T0_1~Sex*Diet*+Litter_Size+Family,sF16Phenotypes))

sF16Phenotypes$T15_1[as.numeric(rownames(data.frame(residuals(aov(T15_1~Sex*Diet*+Litter_Size+Family,sF16Phenotypes)))))]<-residuals(aov(T15_1~Sex*Diet*+Litter_Size+Family,sF16Phenotypes))

sF16Phenotypes$T30_1[as.numeric(rownames(data.frame(residuals(aov(T30_1~Sex*Diet*+Litter_Size+Family,sF16Phenotypes)))))]<-residuals(aov(T30_1~Sex*Diet*+Litter_Size+Family,sF16Phenotypes))

sF16Phenotypes$T60_1[as.numeric(rownames(data.frame(residuals(aov(T60_1~Sex*Diet*+Litter_Size+Family,sF16Phenotypes)))))]<-residuals(aov(T60_1~Sex*Diet*+Litter_Size+Family,sF16Phenotypes))

sF16Phenotypes$T120_1[as.numeric(rownames(data.frame(residuals(aov(T120_1~Sex*Diet*+Litter_Size+Family,sF16Phenotypes)))))]<-residuals(aov(T120_1~Sex*Diet*+Litter_Size+Family,sF16Phenotypes))

sF16Phenotypes$AUC1[as.numeric(rownames(data.frame(residuals(aov(AUC1~Sex*Diet*+Litter_Size+Family,sF16Phenotypes)))))]<-residuals(aov(AUC1~Sex*Diet*+Litter_Size+Family,sF16Phenotypes))


sF16Phenotypes$T0_2[as.numeric(rownames(data.frame(residuals(aov(T0_2~Sex*Diet*+Litter_Size+Family,sF16Phenotypes)))))]<-residuals(aov(T0_2~Sex*Diet*+Litter_Size+Family,sF16Phenotypes))

sF16Phenotypes$T15_2[as.numeric(rownames(data.frame(residuals(aov(T15_2~Sex*Diet*+Litter_Size+Family,sF16Phenotypes)))))]<-residuals(aov(T15_2~Sex*Diet*+Litter_Size+Family,sF16Phenotypes))

sF16Phenotypes$T30_2[as.numeric(rownames(data.frame(residuals(aov(T30_2~Sex*Diet*+Litter_Size+Family,sF16Phenotypes)))))]<-residuals(aov(T30_2~Sex*Diet*+Litter_Size+Family,sF16Phenotypes))

sF16Phenotypes$T60_2[as.numeric(rownames(data.frame(residuals(aov(T60_2~Sex*Diet*+Litter_Size+Family,sF16Phenotypes)))))]<-residuals(aov(T60_2~Sex*Diet*+Litter_Size+Family,sF16Phenotypes))

sF16Phenotypes$T120_2[as.numeric(rownames(data.frame(residuals(aov(T120_2~Sex*Diet*+Litter_Size+Family,sF16Phenotypes)))))]<-residuals(aov(T120_2~Sex*Diet*+Litter_Size+Family,sF16Phenotypes))

sF16Phenotypes$AUC2[as.numeric(rownames(data.frame(residuals(aov(AUC2~Sex*Diet*+Litter_Size+Family,sF16Phenotypes)))))]<-residuals(aov(AUC2~Sex*Diet*+Litter_Size+Family,sF16Phenotypes))


sF16Phenotypes$NecWeight[as.numeric(rownames(data.frame(residuals(aov(NecWeight~Sex*Diet*+Litter_Size+Family,sF16Phenotypes)))))]<-residuals(aov(NecWeight~Sex*Diet*+Litter_Size+Family,sF16Phenotypes))

sF16Phenotypes$Heart[as.numeric(rownames(data.frame(residuals(aov(Heart~Sex*Diet*+Litter_Size+Family,sF16Phenotypes)))))]<-residuals(aov(Heart~Sex*Diet*+Litter_Size+Family,sF16Phenotypes))

sF16Phenotypes$L.Kid[as.numeric(rownames(data.frame(residuals(aov(L.Kid~Sex*Diet*+Litter_Size+Family,sF16Phenotypes)))))]<-residuals(aov(L.Kid~Sex*Diet*+Litter_Size+Family,sF16Phenotypes))

sF16Phenotypes$R.Kid[as.numeric(rownames(data.frame(residuals(aov(R.Kid~Sex*Diet*+Litter_Size+Family,sF16Phenotypes)))))]<-residuals(aov(R.Kid~Sex*Diet*+Litter_Size+Family,sF16Phenotypes))

sF16Phenotypes$Kidney[as.numeric(rownames(data.frame(residuals(aov(Kidney~Sex*Diet*+Litter_Size+Family,sF16Phenotypes)))))]<-residuals(aov(Kidney~Sex*Diet*+Litter_Size+Family,sF16Phenotypes))

sF16Phenotypes$Spleen[as.numeric(rownames(data.frame(residuals(aov(Spleen~Sex*Diet*+Litter_Size+Family,sF16Phenotypes)))))]<-residuals(aov(Spleen~Sex*Diet*+Litter_Size+Family,sF16Phenotypes))

sF16Phenotypes$Liver[as.numeric(rownames(data.frame(residuals(aov(Liver~Sex*Diet*+Litter_Size+Family,sF16Phenotypes)))))]<-residuals(aov(Liver~Sex*Diet*+Litter_Size+Family,sF16Phenotypes))

sF16Phenotypes$Go.dal[as.numeric(rownames(data.frame(residuals(aov(Go.dal~Sex*Diet*+Litter_Size+Family,sF16Phenotypes)))))]<-residuals(aov(Go.dal~Sex*Diet*+Litter_Size+Family,sF16Phenotypes))

sF16Phenotypes$Kid_Fat[as.numeric(rownames(data.frame(residuals(aov(Kid_Fat~Sex*Diet*+Litter_Size+Family,sF16Phenotypes)))))]<-residuals(aov(Kid_Fat~Sex*Diet*+Litter_Size+Family,sF16Phenotypes))

sF16Phenotypes$Mesen[as.numeric(rownames(data.frame(residuals(aov(Mesen~Sex*Diet*+Litter_Size+Family,sF16Phenotypes)))))]<-residuals(aov(Mesen~Sex*Diet*+Litter_Size+Family,sF16Phenotypes))

sF16Phenotypes$Ing[as.numeric(rownames(data.frame(residuals(aov(Ing~Sex*Diet*+Litter_Size+Family,sF16Phenotypes)))))]<-residuals(aov(Ing~Sex*Diet*+Litter_Size+Family,sF16Phenotypes))

sF16Phenotypes$FFA[as.numeric(rownames(data.frame(residuals(aov(FFA~Sex*Diet*+Litter_Size+Family,sF16Phenotypes)))))]<-residuals(aov(FFA~Sex*Diet*+Litter_Size+Family,sF16Phenotypes))

sF16Phenotypes$GLC[as.numeric(rownames(data.frame(residuals(aov(GLC~Sex*Diet*+Litter_Size+Family,sF16Phenotypes)))))]<-residuals(aov(GLC~Sex*Diet*+Litter_Size+Family,sF16Phenotypes))

sF16Phenotypes$CHOL[as.numeric(rownames(data.frame(residuals(aov(CHOL~Sex*Diet*+Litter_Size+Family,sF16Phenotypes)))))]<-residuals(aov(CHOL~Sex*Diet*+Litter_Size+Family,sF16Phenotypes))

sF16Phenotypes$TG[as.numeric(rownames(data.frame(residuals(aov(TG~Sex*Diet*+Litter_Size+Family,sF16Phenotypes)))))]<-residuals(aov(TG~Sex*Diet*+Litter_Size+Family,sF16Phenotypes))

sF16Phenotypes$Insulin[as.numeric(rownames(data.frame(residuals(aov(Insulin~Sex*Diet*+Litter_Size+Family,sF16Phenotypes)))))]<-residuals(aov(Insulin~Sex*Diet*+Litter_Size+Family,sF16Phenotypes))

sF16Phenotypes$Tot_Fat[as.numeric(rownames(data.frame(residuals(aov(Tot_Fat~Sex*Diet*+Litter_Size+Family,sF16Phenotypes)))))]<-residuals(aov(Tot_Fat~Sex*Diet*+Litter_Size+Family,sF16Phenotypes))

unique(apply(sF16Phenotypes[7:8],1,function(S) paste(S,collapse=".")))

}















Bar_Dir<-"ExP_Epistasis_Barplots"
dir.create(Bar_Dir)



# # ### Epistasis Test Functions
# {

# AddAndImpTest<-function(locus1,locus2,ImprintingScores,AdditiveScores,Phenotypes,phenOfInterest,DEgene,ASEgene,RANDOM=FALSE){
	
	# ### Pulling matched phenotypes for individuals from Imprinting scores object
	# matchedPhenotypes<-Phenotypes[match(colnames(ImprintingScores[locus1,8:ncol(ImprintingScores)]),Phenotypes$Animal),]
	
	# ### Pulling matched columns for phenotypes of interest
	# matchedPhenIndex<-match(phenOfInterest,colnames(matchedPhenotypes))
	
	# ### Construct data frame from marker genotypes and matched phenotypes
	# datMat<-cbind(t(data.frame(ImprintingScores[c(locus1,locus2),8:ncol(ImprintingScores)])), matchedPhenotypes[,matchedPhenIndex])
	
	
	# colnames(datMat)<-c("L1_Imp", "L2_Imp","phen")
	# datMat<-data.frame(apply(datMat,2,as.numeric))
	
	# if(RANDOM == TRUE){
		# datMat$phen <-sample(datMat$phen,nrow(datMat))
	# }
	
	
	
	# ### Through out any animals with no phenotype
	# datMat<-datMat[!is.na(datMat$phen),]
	
	# ### Run ANOVA model and return F-value's 
	# c( summary(aov(phen~.*.,data.frame(datMat)))[[1]][["Pr(>F)"]][1:3],DEgene,ASEgene,phenOfInterest,ImprintingScores[locus1,4],ImprintingScores[locus2,4] )
	
# }





# ### This runs the test for All combinations of a gene Pairs markers for the relavent phenotype
# genePairEpiTest<-function(gDE,gImp,SNPS,ImprintingScores,AdditiveScores,Phenotypes,phenOfInterest,RANDOM=FALSE){

	# ### Pulls all loci for DE Gene
	# gDELoci<-match(unique(subset(SNPS,Gene==gDE)$LocusID), ImprintingScores$LocusID)
	
	# ### Pulls all loci for ASE Gene
	# gImpLoci<-match(unique(subset(SNPS,Gene== gImp)$LocusID), ImprintingScores$LocusID)
	
	# ### Generates matrix of All Loci Pairs
	# LociMat<-matrix(unlist(lapply(gDELoci, function(A) lapply(gImpLoci,function(B) return(c(A,B))) )),byrow=TRUE,ncol=2)
	
	# ### Runs Test on all pairs in matrix. Returns list of F-statistics
	# unlist((apply(LociMat,1,function(X) AddAndImpTest(X[1], X[2], ImprintingScores, AdditiveScores, Phenotypes, phenOfInterest, gDE, gImp, RANDOM))))

# }






# PlotgenePairEpiTest<-function(gDE,gImp,SNPS,ImprintingScores,AdditiveScores,Phenotypes,phenOfInterest, BARDIR){

	# ### Pulls all loci for DE Gene
	# gDELoci<-match(unique(subset(SNPS,Gene==gDE)$LocusID), ImprintingScores$LocusID)
	
	# ### Pulls all loci for ASE Gene
	# gImpLoci<-match(unique(subset(SNPS,Gene== gImp)$LocusID), ImprintingScores$LocusID)
	
	# ### Generates matrix of All Loci Pairs
	# LociMat<-matrix(unlist(lapply(gDELoci, function(A) lapply(gImpLoci,function(B) return(c(A,B))) )),byrow=TRUE,ncol=2)
	
	# ### Runs Test on all pairs in matrix. Returns list of F-statistics
	# unlist((apply(LociMat,1,function(X) PlotImpTest(X[1], X[2], ImprintingScores, AdditiveScores, Phenotypes, phenOfInterest, gDE, gImp, BARDIR))))

# }





# }










PlotImpTest<-function(locus1,locus2,ImprintingScores,AdditiveScores,Phenotypes,phenOfInterest,DEgene,ASEgene,BARDIR){
	
	### Pulling matched phenotypes for individuals from Imprinting scores object
	matchedPhenotypes<-Phenotypes[match(colnames(ImprintingScores[locus1,8:ncol(ImprintingScores)]),Phenotypes$Animal),]
	
	### Pulling matched columns for phenotypes of interest
	matchedPhenIndex<-match(phenOfInterest,colnames(matchedPhenotypes))
	
	### Construct data frame from marker genotypes and matched phenotypes
	# datMat<-cbind(t(data.frame(ImprintingScores[c(locus1,locus2),8:ncol(ImprintingScores)])), matchedPhenotypes[,matchedPhenIndex])
	
	datMat<-cbind(
	
	t(data.frame(ImprintingScores[c(locus1,locus2),8:ncol(ImprintingScores)])),
	t(data.frame(AdditiveScores[c(locus1,locus2),8:ncol(AdditiveScores)])),
	matchedPhenotypes[,matchedPhenIndex]
	
	)

	
	print(ImprintingScores[c(locus1,locus2),4])
	
	
	colnames(datMat)<-c("L1_Imp", "L2_Imp","L1_Add", "L2_Add","phen")
	datMat<-data.frame(apply(datMat,2,as.numeric))
	
	# LL i:0 a:1
	# LS i:1 a:0
	# SL i:-1 a:0
	# SS i:0 a:-1
	
	datMat$L1_Genotype<-NA
	
	datMat$L1_Genotype[datMat$L1_Imp==0 & datMat$L1_Add==1]<-"LL"
	datMat$L1_Genotype[datMat$L1_Imp==1 & datMat$L1_Add==0]<-"LS"
	datMat$L1_Genotype[datMat$L1_Imp==-1 & datMat$L1_Add==0]<-"SL"
	datMat$L1_Genotype[datMat$L1_Imp==0 & datMat$L1_Add==-1]<-"SS"
	
	datMat$L1_PaternalAllele<-NA
	datMat$L1_PaternalAllele[datMat$L1_Genotype %in% c("LL","SL") ]<-"LG/J"
	datMat$L1_PaternalAllele[datMat$L1_Genotype %in% c("SS","LS") ]<-"SM/J"
	
	
	datMat$L2_Genotype<-NA
		
	datMat$L2_Genotype[datMat$L2_Imp==0 & datMat$L2_Add==1]<-"LL"
	datMat$L2_Genotype[datMat$L2_Imp==1 & datMat$L2_Add==0]<-"LS"
	datMat$L2_Genotype[datMat$L2_Imp==-1 & datMat$L2_Add==0]<-"SL"
	datMat$L2_Genotype[datMat$L2_Imp==0 & datMat$L2_Add==-1]<-"SS"
	
	
	datMat$L2_PaternalAllele<-NA
	datMat$L2_PaternalAllele[datMat$L2_Genotype %in% c("LL","SL") ]<-"LG/J"
	datMat$L2_PaternalAllele[datMat$L2_Genotype %in% c("SS","LS") ]<-"SM/J"

	
	datMat<-subset(datMat, !is.na(L1_Genotype) & !is.na(L2_Genotype))
	head(datMat)
	
	datInt <- ddply(subset(datMat, L1_Imp %in% c(-1,1) & L2_Imp %in% c(-1,1) ), .(L1_Genotype, L2_Genotype), summarise, val = median(phen, na.rm=TRUE))
	
	print(head(datMat))
	
	print(datInt)
	 
	PLOT<-ggplot(subset(datMat, L1_Imp %in% c(-1,1) & L2_Imp %in% c(-1,1) ),aes(x=as.factor(L1_Genotype),y=phen,fill=as.factor(L2_Genotype)))+
	geom_boxplot()+
	geom_line( data = datInt, aes(y = val, group=as.factor(L2_Genotype)),lwd=1,position=position_dodge(width=0.75) )
	
	# PLOT<-ggplot(subset(datMat, L1_Imp %in% c(-1,1) & L2_Imp %in% c(-1,1) ),aes(x=as.factor(L1_Genotype),y=phen,fill=as.factor(L2_Genotype)))+
	# geom_boxplot()+
	# # geom_jitter(size=0.3,width=0.2)+
	# facet_grid(rows=vars(L2_Genotype))
	
	# PLOT<-ggplot(datMat,aes(x=as.factor(L1_PaternalAllele),y=phen))+
	# geom_boxplot()+
	# geom_jitter(size=0.3,width=0.2)+
	# facet_grid(cols=vars(L2_PaternalAllele))
	print(datInt)

	return(PLOT)	
	# print(PLOT)
	
	### Through out any animals with no phenotype
	datMat<-datMat[!is.na(datMat$phen),]
	
	
	DFMEANS<-data.frame(
	L1_Imp=c(1,-1,1,-1),	
	L2_Imp=c(1,1,-1,-1),
	Means=c(
	mean(subset(datMat, L1_Imp == 1 & L2_Imp == 1 )$phen),
	mean(subset(datMat, L1_Imp == -1 & L2_Imp == 1 )$phen),
	mean(subset(datMat, L1_Imp == 1 & L2_Imp == -1 )$phen),
	mean(subset(datMat, L1_Imp == -1 & L2_Imp == -1 )$phen)
	)
	
	)

	
	# print(datMat)
	

	
	
	
	# FILENAME<-paste(BARDIR,"/EpistatsisPlot","_",phenOfInterest,"_",DEgene,"_",ASEgene,"_",ImprintingScores$Marker.ID[locus1],"_",ImprintingScores$Marker.ID[locus2],".tiff", sep="")
	
	
	# PLOT<-ggplot(DFMEANS,aes(x=L1_Imp, y=Means))+
	# facet_grid(rows = vars(L2_Imp) )+
	# geom_bar(stat='identity')+
	# scale_x_discrete()+
	# ggtitle( paste(phenOfInterest," ~ ",DEgene," X ", ASEgene,"\n",ImprintingScores$Marker.ID[locus1]," : ",ImprintingScores$Marker.ID[locus2], sep="") )+
	# xlab("DE Genotype")+
	# ylab("Phenotype")
	
	# ggsave(as.character(FILENAME), PLOT, units="in",width=3.25,height=3)
	
	### Run ANOVA model and return F-value's 
	# c( summary(aov(phen~.*.,data.frame(datMat)))[[1]][["Pr(>F)"]][1:3],DEgene,ASEgene,phenOfInterest,ImprintingScores[locus1,4],ImprintingScores[locus2,4] )
	
}





PHENOTYPES<-sF16Phenotypes
PHENOTYPES<-F16Phenotypes

IMPRINTINGSCORES<-NetImprintingScores
ADDITIVESCORES<-NetAdditiveScores
MEDIANS<-NetMedian

### Subset network object to DE~ASE~Phen sets where both the ASE and DE gene has associated markers
Net<-Net[(Net$ASE.Gene %in% as.character(unique(GeneSNPs$Gene)) & Net$Target.Gene %in% as.character(unique(GeneSNPs$Gene))),]












#Nnat
(1:nrow(IMPRINTINGSCORES))[IMPRINTINGSCORES$LocusID=="2_157295559"]
(1:nrow(IMPRINTINGSCORES))[IMPRINTINGSCORES$Marker.ID=="UT-2-158.095429"]



#F2r
(1:nrow(IMPRINTINGSCORES))[IMPRINTINGSCORES$LocusID=="13_98551175"]
(1:nrow(IMPRINTINGSCORES))[IMPRINTINGSCORES$Marker.ID=="wu-rs13481958"]

#Cdkn1c
(1:nrow(IMPRINTINGSCORES))[IMPRINTINGSCORES$LocusID=="7_150663801"]


#Plcd1
(1:nrow(IMPRINTINGSCORES))[IMPRINTINGSCORES$LocusID=="9_117827881"]




Plot.Nnat.F2r.T0_2<-PlotImpTest(34,15, IMPRINTINGSCORES, ADDITIVESCORES, PHENOTYPES,as.character("T0_2"),"F2r","Nnat", Bar_Dir)


Plot.Nnat.F2r.T0_2+
ylab("Basal Glucose")+
xlab("F2r Genotype")+
labs(fill = "Nnat\nGenotype")+
theme(
	panel.background = element_rect(fill="white"),
	panel.grid.major=element_line(color="gray90", linetype="solid"), #adds major grid lines
	panel.border=element_rect(fill=NA, color="black", size=1.5, linetype="solid"), #draws a black border around plot
	axis.line=element_blank(),
	)+
theme(legend.title=element_text(size=8))





PlotImpTest(35,21, IMPRINTINGSCORES, ADDITIVESCORES, PHENOTYPES,as.character("T0_2"),"F2r","Cdkn1c", Bar_Dir)



PlotImpTest(35,28, IMPRINTINGSCORES, ADDITIVESCORES, PHENOTYPES,as.character("T0_2"),"F2r","Plcd1", Bar_Dir)




