

rm(list = ls())
set.seed(0)
require("edgeR")
library(lsr)

library(RUVSeq)

library(ggplot2)
library(cowplot)
library(MASS)

### Functions



IDCon<-read.delim("IDconversions_LGSM.txt")

IDCon$Gene.name<-as.character(IDCon$Gene.name)



### Import Data

#RawData<-read.delim("All_WAT_Normed.txt",header=TRUE)

RawData<-read.delim("Normalized_Counts.tsv",header=TRUE)
Data<-RawData[,2:ncol(RawData)]
rownames(Data)<-as.character(RawData$Ensembl)


SampleInfo<-read.delim("All_BAT_WAT_SampleInfo.txt",header=TRUE)
SampleInfo<-subset(SampleInfo,Tissue=="r")

#Order Sample Info
SampleInfo<-SampleInfo[match(colnames(Data),as.character(SampleInfo$GTAC.Tag)),]


Sex<-SampleInfo$Sex
Diet<-SampleInfo$Diet
Age<-SampleInfo$Age
Genotype<-SampleInfo$Genotype


#Sex<-as.numeric(Sex)-1
#Age<-as.numeric(as.factor(Age))-1
#Diet<-as.numeric(Diet)-1
Counts.Mat<-matrix(unlist(Data),ncol=ncol(Data))
colnames(Counts.Mat)<-colnames(Data)
rownames(Counts.Mat)<-rownames(Data)


colnames(Counts.Mat)==as.character(SampleInfo$GTAC.Tag)

#Counts.Mat<-Counts.Mat[,Age==20] ### Remove 30 week animals

predictors<-data.frame(Diet,Sex,Age,Genotype)

dim(Data)


Exclude<-NULL
Exclude<-c("17","39")

subset(SampleInfo,Age==20)[as.numeric(Exclude),]

### Nnat
NNAT<-data.frame(Counts.Mat[match("ENSMUSG00000067786",rownames(Counts.Mat)),],predictors)

colnames(NNAT)<-c("Count","DIET","SEX","AGE","GENOTYPE")

NNAT<-subset(NNAT,AGE==20)

NNAT<-NNAT[!(1:nrow(NNAT) %in% Exclude),]



#NNAT$GENOTYPE[NNAT$GENOTYPE=="LS"]<-"SS"
#NNAT$GENOTYPE[NNAT$GENOTYPE=="SL"]<-"LL"
NNAT$GENOTYPE<-factor(as.character(NNAT$GENOTYPE))



summary(aov(Count~GENOTYPE*AGE, NNAT))
summary(aov(Count~GENOTYPE*AGE*DIET, NNAT))
summary(aov(Count~GENOTYPE*AGE*DIET*SEX, NNAT))

tiff(paste("Nnat_Genotype.tiff",sep=""),units='in',width=5.5,height=5.5,res=500)

ggplot(NNAT, aes(x=GENOTYPE, y=Count),notch=TRUE) + 
geom_boxplot(fill='#A4A4A4',outlier.colour ="white") +
labs(x="Genotype",y="Normalized Expression") + 
theme(
    axis.text.x=element_text(angle=0, size=15),
    axis.text.y=element_text(angle=0, size=15),
    axis.title.x=element_text(angle=0, face='bold', size=30),
    axis.title.y=element_text(angle=90, face='bold', size=30)
)

dev.off()


tiff(paste("Nnat_Genotype_Diet.tiff",sep=""),units='in',width=5.5,height=5.5,res=500)

ggplot(NNAT, aes(x=GENOTYPE, y=Count),notch=TRUE) + 
geom_boxplot(fill='#A4A4A4',outlier.colour ="white") +
labs(x="Genotype",y="Normalized Expression") + 
theme(
    axis.text.x=element_text(angle=0, size=15),
    axis.text.y=element_text(angle=0, size=15),
    axis.title.x=element_text(angle=0, face='bold', size=30),
    axis.title.y=element_text(angle=90, face='bold', size=30)
) + facet_wrap(~DIET)

dev.off()

tiff(paste("Nnat_Genotype_Diet_Sex.tiff",sep=""),units='in',width=5.5,height=5.5,res=500)

ggplot(NNAT, aes(x=GENOTYPE, y=Count),notch=TRUE) + 
geom_boxplot(fill='#A4A4A4',outlier.colour ="white") +
labs(x="Genotype",y="Normalized Expression") + 
theme(
    axis.text.x=element_text(angle=0, size=15),
    axis.text.y=element_text(angle=0, size=15),
    axis.title.x=element_text(angle=0, face='bold', size=30),
    axis.title.y=element_text(angle=90, face='bold', size=30)
) + facet_wrap(~DIET*SEX)

dev.off()

F.NNAT<-subset(NNAT,SEX=="Female")
HF.NNAT <-subset(NNAT,SEX=="Female" & DIET=="Highfat")

minY<-min(NNAT$Count)*0.9
maxY<-max(NNAT$Count)*1.1

breakSizeY<-round(round((maxY-minY)/8)/5)*5

tiff(paste("Nnat_Genotype_Diet_Females.tiff",sep=""),units='in',width=5.5,height=5.5,res=500)

ggplot(F.NNAT, aes(x= GENOTYPE, y=Count),notch=TRUE) + 
geom_boxplot(fill='#A4A4A4',outlier.colour ="white") +
labs(x="Genotype",y="Normalized Expression") + 
theme(
    axis.text.x=element_text(angle=0, size=15),
    axis.text.y=element_text(angle=0, size=15),
    axis.title.x=element_text(angle=0, face='bold', size=30),
    axis.title.y=element_text(angle=90, face='bold', size=30)
) + facet_wrap(~DIET)

dev.off()


tiff(paste("Nnat_Genotype_HighFat_Females.tiff",sep=""),units='in',width=5.5,height=5.5,res=500)

ggplot(HF.NNAT, aes(x= GENOTYPE, y=Count),notch=TRUE) + 
geom_boxplot(fill='#A4A4A4',outlier.colour ="white") +
labs(x="Genotype",y="Normalized Expression") + 
theme(
    axis.text.x=element_text(angle=0, size=15),
    axis.text.y=element_text(angle=0, size=15),
    axis.title.x=element_text(angle=0, face='bold', size=30),
    axis.title.y=element_text(angle=90, face='bold', size=30)
) + facet_wrap(~DIET)

dev.off()










### ATF4
ATF4<-data.frame(Counts.Mat[match("ENSMUSG00000042406",rownames(Counts.Mat)),],predictors)

colnames(ATF4)<-c("Count","DIET","SEX","AGE","GENOTYPE")

ATF4<-subset(ATF4,AGE==20)

ATF4<-ATF4[!(1:nrow(ATF4) %in% Exclude),]



#ATF4$GENOTYPE[ATF4$GENOTYPE=="LS"]<-"SS"
#ATF4$GENOTYPE[ATF4$GENOTYPE=="SL"]<-"LL"
ATF4$GENOTYPE<-factor(as.character(ATF4$GENOTYPE))



summary(aov(Count~GENOTYPE*AGE, ATF4))
summary(aov(Count~GENOTYPE*AGE*DIET, ATF4))
summary(aov(Count~GENOTYPE*AGE*DIET*SEX, ATF4))

tiff(paste("ATF4_Genotype.tiff",sep=""),units='in',width=5.5,height=5.5,res=500)

ggplot(ATF4, aes(x=GENOTYPE, y=Count),notch=TRUE) + 
geom_boxplot(fill='#A4A4A4',outlier.colour ="white") +
labs(x="Genotype",y="Normalized Expression") + 
theme(
    axis.text.x=element_text(angle=0, size=15),
    axis.text.y=element_text(angle=0, size=15),
    axis.title.x=element_text(angle=0, face='bold', size=30),
    axis.title.y=element_text(angle=90, face='bold', size=30)
)

dev.off()


tiff(paste("ATF4_Genotype_Diet.tiff",sep=""),units='in',width=5.5,height=5.5,res=500)

ggplot(ATF4, aes(x=GENOTYPE, y=Count),notch=TRUE) + 
geom_boxplot(fill='#A4A4A4',outlier.colour ="white") +
labs(x="Genotype",y="Normalized Expression") + 
theme(
    axis.text.x=element_text(angle=0, size=15),
    axis.text.y=element_text(angle=0, size=15),
    axis.title.x=element_text(angle=0, face='bold', size=30),
    axis.title.y=element_text(angle=90, face='bold', size=30)
) + facet_wrap(~DIET)

dev.off()

tiff(paste("ATF4_Genotype_Diet_Sex.tiff",sep=""),units='in',width=5.5,height=5.5,res=500)

ggplot(ATF4, aes(x=GENOTYPE, y=Count),notch=TRUE) + 
geom_boxplot(fill='#A4A4A4',outlier.colour ="white") +
labs(x="Genotype",y="Normalized Expression") + 
theme(
    axis.text.x=element_text(angle=0, size=15),
    axis.text.y=element_text(angle=0, size=15),
    axis.title.x=element_text(angle=0, face='bold', size=30),
    axis.title.y=element_text(angle=90, face='bold', size=30)
) + facet_wrap(~DIET*SEX)

dev.off()

F.ATF4<-subset(ATF4,SEX=="Female")
HF.ATF4 <-subset(ATF4,SEX=="Female" & DIET=="Highfat")

minY<-min(ATF4$Count)*0.9
maxY<-max(ATF4$Count)*1.1

breakSizeY<-round(round((maxY-minY)/8)/5)*5

tiff(paste("ATF4_Genotype_Diet_Females.tiff",sep=""),units='in',width=5.5,height=5.5,res=500)

ggplot(F.ATF4, aes(x= GENOTYPE, y=Count),notch=TRUE) + 
geom_boxplot(fill='#A4A4A4',outlier.colour ="white") +
labs(x="Genotype",y="Normalized Expression") + 
theme(
    axis.text.x=element_text(angle=0, size=15),
    axis.text.y=element_text(angle=0, size=15),
    axis.title.x=element_text(angle=0, face='bold', size=30),
    axis.title.y=element_text(angle=90, face='bold', size=30)
) + facet_wrap(~DIET)

dev.off()


tiff(paste("ATF4_Genotype_HighFat_Females.tiff",sep=""),units='in',width=5.5,height=5.5,res=500)

ggplot(HF.ATF4, aes(x= GENOTYPE, y=Count),notch=TRUE) + 
geom_boxplot(fill='#A4A4A4',outlier.colour ="white") +
labs(x="Genotype",y="Normalized Expression") + 
theme(
    axis.text.x=element_text(angle=0, size=15),
    axis.text.y=element_text(angle=0, size=15),
    axis.title.x=element_text(angle=0, face='bold', size=30),
    axis.title.y=element_text(angle=90, face='bold', size=30)
) + facet_wrap(~DIET)

dev.off()













### Mogat1

MOGAT1<-data.frame(Counts.Mat[match("ENSMUSG00000012187",rownames(Counts.Mat)),],predictors)

colnames(MOGAT1)<-c("Count","DIET","SEX","AGE","GENOTYPE")

MOGAT1<-subset(MOGAT1,AGE==20)

MOGAT1<-MOGAT1[!(1:nrow(MOGAT1) %in% Exclude),]


#MOGAT1$GENOTYPE[MOGAT1$GENOTYPE=="LS"]<-"SS"
#MOGAT1$GENOTYPE[MOGAT1$GENOTYPE=="SL"]<-"LL"
MOGAT1$GENOTYPE<-factor(as.character(MOGAT1$GENOTYPE))



summary(aov(Count~GENOTYPE*AGE, MOGAT1))
summary(aov(Count~GENOTYPE*AGE*DIET, MOGAT1))
summary(aov(Count~GENOTYPE*AGE*DIET*SEX, MOGAT1))


tiff(paste("Mogat1_Genotype.tiff",sep=""),units='in',width=5.5,height=5.5,res=500)

ggplot(MOGAT1, aes(x=GENOTYPE, y=Count),notch=TRUE) + 
geom_boxplot(fill='#A4A4A4',outlier.colour ="white") +
labs(x="Genotype",y="Normalized Expression") + 
theme(
    axis.text.x=element_text(angle=0, size=15),
    axis.text.y=element_text(angle=0, size=15),
    axis.title.x=element_text(angle=0, face='bold', size=30),
    axis.title.y=element_text(angle=90, face='bold', size=30)
)

dev.off()

tiff(paste("Mogat1_Genotype_Diet.tiff",sep=""),units='in',width=5.5,height=5.5,res=500)

ggplot(MOGAT1, aes(x=GENOTYPE, y=Count),notch=TRUE) + 
geom_boxplot(fill='#A4A4A4',outlier.colour ="white") +
labs(x="Genotype",y="Normalized Expression") + 
theme(
    axis.text.x=element_text(angle=0, size=15),
    axis.text.y=element_text(angle=0, size=15),
    axis.title.x=element_text(angle=0, face='bold', size=30),
    axis.title.y=element_text(angle=90, face='bold', size=30)
) + facet_wrap(~DIET)

dev.off()



tiff(paste("Mogat1_Genotype_Diet_Sex.tiff",sep=""),units='in',width=5.5,height=5.5,res=500)

ggplot(MOGAT1, aes(x=GENOTYPE, y=Count),notch=TRUE) + 
geom_boxplot(fill='#A4A4A4',outlier.colour ="white") +
labs(x="Genotype",y="Normalized Expression") + 
theme(
    axis.text.x=element_text(angle=0, size=15),
    axis.text.y=element_text(angle=0, size=15),
    axis.title.x=element_text(angle=0, face='bold', size=30),
    axis.title.y=element_text(angle=90, face='bold', size=30)
) + facet_wrap(~DIET*SEX)

dev.off()

F.MOGAT1<-subset(MOGAT1,SEX=="Female")
HF.MOGAT1<-subset(MOGAT1,SEX=="Female" & DIET=="Highfat")

minY<-min(MOGAT1$Count)*0.9
maxY<-max(MOGAT1$Count)*1.1

breakSizeY<-round(round((maxY-minY)/8)/5)*5


tiff(paste("Mogat1_Genotype_Diet_Females.tiff",sep=""),units='in',width=5.5,height=5.5,res=500)

ggplot(F.MOGAT1, aes(x= GENOTYPE, y=Count),notch=TRUE) + 
geom_boxplot(fill='#A4A4A4',outlier.colour ="white") +
labs(x="Genotype",y="Normalized Expression") + 
theme(
    axis.text.x=element_text(angle=0, size=15),
    axis.text.y=element_text(angle=0, size=15),
    axis.title.x=element_text(angle=0, face='bold', size=30),
    axis.title.y=element_text(angle=90, face='bold', size=30)
) + facet_wrap(~DIET)

dev.off()


tiff(paste("Mogat1_Genotype_HighFat_Females.tiff",sep=""),units='in',width=5.5,height=5.5,res=500)

ggplot(HF.MOGAT1, aes(x= GENOTYPE, y=Count),notch=TRUE) + 
geom_boxplot(fill='#A4A4A4',outlier.colour ="white") +
labs(x="Genotype",y="Normalized Expression") + 
theme(
    axis.text.x=element_text(angle=0, size=15),
    axis.text.y=element_text(angle=0, size=15),
    axis.title.x=element_text(angle=0, face='bold', size=30),
    axis.title.y=element_text(angle=90, face='bold', size=30)
) + facet_wrap(~DIET)

dev.off()









DataFrame<-data.frame(MOGAT1$Count,NNAT$Count,MOGAT1[,2:5])
colnames(DataFrame)<-c("Mogat1","Nnat","DIET","SEX","AGE","GENOTYPE")

#DataFrame<-subset(DataFrame,SEX=="Female")
#DataFrame<-subset(DataFrame,SEX=="Female" & DIET=="Highfat")
#DataFrame<-subset(DataFrame, (GENOTYPE=="LS" |GENOTYPE=="SL") )

HF<-subset(DataFrame,DIET=="Highfat" & SEX=="Female")
HM<-subset(DataFrame,DIET=="Highfat" & SEX=="Male")
LF<-subset(DataFrame,DIET=="Lowfat" & SEX=="Female")
LM<-subset(DataFrame,DIET=="Lowfat" & SEX=="Male")

HF.LL<-subset(DataFrame, GENOTYPE=="LL" & DIET=="Highfat" & SEX=="Female")
HF.LS<-subset(DataFrame, GENOTYPE=="LS" & DIET=="Highfat" & SEX=="Female")
HF.SL<-subset(DataFrame, GENOTYPE=="SL" & DIET=="Highfat" & SEX=="Female")
HF.SS<-subset(DataFrame, GENOTYPE=="SS" & DIET=="Highfat" & SEX=="Female")

r<-round(cor(DataFrame$Mogat1, DataFrame$Nnat),2)

Rsqd<-c(
round(unlist(summary(lm(HF$Mogat1~HF$Nnat))[9]),2),
round(unlist(summary(lm(HM$Mogat1~HM$Nnat))[9]),2),
round(unlist(summary(lm(LF$Mogat1~LF$Nnat))[9]),2),
round(unlist(summary(lm(LM$Mogat1~LM$Nnat))[9]),2)
)


RsqdLabs<-data.frame(
xp=rep(quantile(DataFrame$Mogat1,0.92),4),
yp=rep(quantile(DataFrame$Nnat,0.05),4),
lab=Rsqd,
DIET=c("Highfat","Highfat","Lowfat","Lowfat"),
SEX=c("Female","Male","Female","Male")
)

cor.test(DataFrame$Mogat1, DataFrame$Nnat)


tiff(paste("Nnat_Mogat_Diet_Sex_Correlation.tiff",sep=""),units='in',width=5.5,height=5.5,res=500)

ggplot(DataFrame, aes(x=Mogat1, y=Nnat)) + 
geom_point(size=1) + 
geom_smooth(method=lm,se=FALSE,size=0.75) + 
labs(x=paste("Mogat1","expression"),y=paste("Nnat","expression")) + 
theme(
    axis.text.x=element_text(angle=0, size=15),
    axis.text.y=element_text(angle=0, size=15),
    axis.title.x=element_text(angle=0, face='bold', size=30),
    axis.title.y=element_text(angle=90, face='bold', size=30)
) + 
facet_wrap(~DIET*SEX) + 
geom_text(aes(
		x=xp,
		y=yp,
		label=paste0('R^2 == ', lab),
		family = "Arial",
		size=10),
		parse = TRUE,
		data=RsqdLabs
)

dev.off()









DataFrame<-data.frame(MOGAT1$Count,ATF4$Count,MOGAT1[,2:5])
colnames(DataFrame)<-c("Mogat1","ATF4","DIET","SEX","AGE","GENOTYPE")

#DataFrame<-subset(DataFrame,SEX=="Female")
#DataFrame<-subset(DataFrame,SEX=="Female" & DIET=="Highfat")
#DataFrame<-subset(DataFrame, (GENOTYPE=="LS" |GENOTYPE=="SL") )

HF<-subset(DataFrame,DIET=="Highfat" & SEX=="Female")
HM<-subset(DataFrame,DIET=="Highfat" & SEX=="Male")
LF<-subset(DataFrame,DIET=="Lowfat" & SEX=="Female")
LM<-subset(DataFrame,DIET=="Lowfat" & SEX=="Male")

r<-round(cor(DataFrame$Mogat1, DataFrame$ATF4),2)

Rsqd<-c(
round(unlist(summary(lm(HF$Mogat1~HF$ATF4))[9]),2),
round(unlist(summary(lm(HM$Mogat1~HM$ATF4))[9]),2),
round(unlist(summary(lm(LF$Mogat1~LF$ATF4))[9]),2),
round(unlist(summary(lm(LM$Mogat1~LM$ATF4))[9]),2)
)


RsqdLabs<-data.frame(
xp=rep(quantile(DataFrame$Mogat1,0.92),4),
yp=rep(quantile(DataFrame$ATF4,0.05),4),
lab=Rsqd,
DIET=c("Highfat","Highfat","Lowfat","Lowfat"),
SEX=c("Female","Male","Female","Male")
)

cor.test(DataFrame$Mogat1, DataFrame$ATF4)


tiff(paste("ATF4_Mogat_Diet_Sex_Correlation.tiff",sep=""),units='in',width=5.5,height=5.5,res=500)

ggplot(DataFrame, aes(x=Mogat1, y=ATF4)) + 
geom_point(size=1) + 
geom_smooth(method=lm,se=FALSE,size=0.75) + 
labs(x=paste("Mogat1","expression"),y=paste("ATF4","expression")) + 
theme(
    axis.text.x=element_text(angle=0, size=15),
    axis.text.y=element_text(angle=0, size=15),
    axis.title.x=element_text(angle=0, face='bold', size=30),
    axis.title.y=element_text(angle=90, face='bold', size=30)
) + 
facet_wrap(~DIET*SEX) + 
geom_text(aes(
		x=xp,
		y=yp,
		label=paste0('R^2 == ', lab),
		family = "Arial",
		size=10),
		parse = TRUE,
		data=RsqdLabs
)

dev.off()






HF.LL<-subset(DataFrame, GENOTYPE=="LL" & DIET=="Highfat" & SEX=="Female")
HF.LS<-subset(DataFrame, GENOTYPE=="LS" & DIET=="Highfat" & SEX=="Female")
HF.SL<-subset(DataFrame, GENOTYPE=="SL" & DIET=="Highfat" & SEX=="Female")
HF.SS<-subset(DataFrame, GENOTYPE=="SS" & DIET=="Highfat" & SEX=="Female")


r<-round(cor(HF$Mogat1, HF$ATF4),2)

Rsqd<-c(
round(unlist(summary(lm(HF.LL$Mogat1~ HF.LL$ATF4))[9]),2),
round(unlist(summary(lm(HF.LS$Mogat1~ HF.LS$ATF4))[9]),2),
round(unlist(summary(lm(HF.SL$Mogat1~ HF.SL$ATF4))[9]),2),
round(unlist(summary(lm(HF.SS$Mogat1~ HF.SS$ATF4))[9]),2)
)


RsqdLabs<-data.frame(
xp=rep(quantile(HF$Mogat1,0.92),4),
yp=rep(quantile(HF$ATF4,0.05),4),
lab=Rsqd,
DIET=c("Highfat","Highfat","Highfat","Highfat"),
SEX=c("Female","Female","Female","Female"),
GENOTYPE=c("LL","LS","SL","SS")
)

cor.test(HF$Mogat1, HF$ATF4)


tiff(paste("ATF4_Mogat_GENOTYPE_HF_Correlation.tiff",sep=""),units='in',width=5.5,height=5.5,res=500)

ggplot(HF, aes(x=Mogat1, y=ATF4)) + 
geom_point(size=1) + 
geom_smooth(method=lm,se=FALSE,size=0.75) + 
labs(x=paste("Mogat1","expression"),y=paste("ATF4","expression")) + 
theme(
    axis.text.x=element_text(angle=0, size=15),
    axis.text.y=element_text(angle=0, size=15),
    axis.title.x=element_text(angle=0, face='bold', size=30),
    axis.title.y=element_text(angle=90, face='bold', size=30)
) + 
facet_wrap(~GENOTYPE) + 
geom_text(aes(
		x=xp,
		y=yp,
		label=paste0('R^2 == ', lab),
		family = "Arial",
		size=10),
		parse = TRUE,
		data=RsqdLabs
)

dev.off()




