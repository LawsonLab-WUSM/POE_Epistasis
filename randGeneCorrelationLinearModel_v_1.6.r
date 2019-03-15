print("Importing Data...")

# [ASE Based Data] [General Data] [ASE File Start] [ASE File stop] [Gen File Start] [Gen File Stop] [Iteration]
args <- commandArgs(trailingOnly = TRUE)

library(MASS)
library(glmnet)
library(methods)
library(boot)


### POE
Data<-read.table(args[1],header=FALSE)
colnames(Data)<-c("Gene","Cross","Sex","Diet","Number","BC","Lbiase","L","S")
Data$Sbiase<-(1-Data$Lbiase)
Data<-subset(Data,Lbiase != "NA") ## Remove NA samples
Data$Lbiase<-Data$L/(Data$L+Data$S)
#Data$Sbiase<-1-Data$Lbiase
Data$Total<-Data$L+Data$S


Data$Cross<-sample(Data$Cross,replace=FALSE)
Data$Sex<-sample(Data$Sex,replace=FALSE)
Data$Diet<-sample(Data$Diet,replace=FALSE)
Data$Total<-sample(Data$Total,replace=FALSE)
Data$Lbiase<-sample(Data$Lbiase,replace=FALSE)


### All
General<-read.table(args[2],header=FALSE)
colnames(General)<-c("Gene","Cross","Sex","Diet","Number","BC","Lbiase","L","S")
General<-subset(General,Lbiase != "NA") ## Remove NA samples
General$Lbiase<-General$L/(General$L+General$S)
#General$Sbiase<-1-General$Lbiase
General$Total<-General$L+General$S


General$Cross<-sample(General$Cross,replace=FALSE)
General$Sex<-sample(General$Sex,replace=FALSE)
General$Diet<-sample(General$Diet,replace=FALSE)
General$Total<-sample(General$Total,replace=FALSE)
General$Lbiase<-sample(General$Lbiase,replace=FALSE)


General<-subset(General, Gene %in% unique(General$Gene)[args[5]:args[6]])
General<-subset(General,Gene %in% rownames(table(General$Gene))[table(General$Gene)>=16])

All_Pass<-scan(file="Passed.txt",what="")
All_Pass<-All_Pass[args[3]:args[4]]

Data<-subset(Data, Gene %in% All_Pass) #Reduce size and search space of Data object

print("Complete")



runLasso2<-function(setA,setB,permNumber){

setB$BT<-setB$Total*((setB$Lbiase-0.5)*2)
setB$MH<-(as.numeric(as.factor(setB$Sex == "M" & setB$Diet == "H"))-1)# 0 or 1
setB$ML<-(as.numeric(as.factor(setB$Sex == "M" & setB$Diet == "L"))-1)
setB$FH<-(as.numeric(as.factor(setB$Sex == "F" & setB$Diet == "H"))-1)
setB$FL<-(as.numeric(as.factor(setB$Sex == "F" & setB$Diet == "L"))-1)
setB$Fd<-as.numeric(as.factor((paste(setB$FH,setB$FL))))-1
setB$Md<-as.numeric(as.factor((paste(setB$MH,setB$ML))))-1
setB$Context<-(as.factor(paste(setB$Sex,setB$Diet,sep="")))


rownames(setB)<-NULL

setB$Total<-scale(setB$Total)
setA$Total<-scale(setA$Total)
setB$BT<-scale(setB$BT)

dat<-data.frame(setB$Lbiase,setB$Total,setB$BT,setB$MH,setB$ML,setB$FH,setB$FL,setB$Fd,setB$Md,setA$Total)
colnames(dat)<-c("Bias","Total","BT","MH","ML","FH","FL","Fd","Md","Y")


set.seed(1)
cv.err<-NULL
Coefs<-NULL
for(z in 1:permNumber){
        train<-sample(nrow(dat),nrow(dat)*0.6,replace = FALSE)
        playMod<-NULL
        #playMod<-glm(Y~Total+BT+Total:MH+Total:ML+Total:FH+Total:FL+BT:MH+BT:ML+BT:FH+BT:FL,dat,family=gaussian)
	playMod<-glm(Y~Total+BT+Total:Fd+Total:Md+BT:Fd+BT:Md,dat,family=gaussian)
        Coefs<-rbind(Coefs,playMod$coefficients)
        cv.err<-c(cv.err,cv.glm(dat,playMod,K=10)$delta[1])
}

mCVE<-round(mean(cv.err),2)
vCVE<-round(sd(cv.err),2)
mCoefs<-apply(Coefs,2,mean)


ACVs<-c(mCoefs,mCVE,vCVE)

return(ACVs)
}

print("Starting Analysis...")
start<-proc.time()
g1<-NULL
g2<-NULL

Tp<-NULL
Bp<-NULL
TBp<-NULL

Te<-NULL
Be<-NULL
TBe<-NULL

minN<-16
l<-1
Pip<-NULL
Pfp<-NULL
k<-1


AllVals<-NULL
NumberofIterationsLASSO<-10


##### dont forget to change back to general!
for(gA in unique(General$Gene)){ ### All_Pass is from Interactions script
	gAdat<-subset(General,Gene==gA)
	Pip[k]<-l
	for(gB in All_Pass){
		hold<-NULL
		if(gA != gB) gBdat<-subset(Data,Gene==gB) else gBdat<-data.frame(1,1)
		
		if(gA != gB & minN<=nrow(gAdat) & minN<=nrow(gBdat)) setA<-gAdat[match(intersect(gBdat$Number,gAdat$Number),gAdat$Number),] else setA<-data.frame(1,1)
		if(gA != gB & minN<=nrow(gAdat) & minN<=nrow(gBdat)) setB<-gBdat[match(intersect(gAdat$Number,gBdat$Number),gBdat$Number),] else setB<-data.frame(1,1)
		gBdat<-NULL
		if(gA != gB) pass<-(minN<=nrow(setA) & minN<=nrow(setB)) else pass<-FALSE

		if(pass) hold<-runLasso2(setA,setB,NumberofIterationsLASSO)

		#if(pass) predY<-hold[11]+hold[12]*as.numeric(setB$Sex)+hold[13]*as.numeric(setB$Diet)+hold[14]*setB$Lbiase+hold[15]*setB$Total
		#if(pass) yBar<-mean(setA$Total)
		#if(pass) SST<-sum((setA$Total-yBar)^2)
		#if(pass) SSR<-sum((setA$Total-predY)^2)
		#if(pass) R2<- 1-(SSR/SST) else R2<-0
		
		if(pass) g1[l]<-gA
		if(pass) g2[l]<-gB
		if(pass) AllVals<-rbind(AllVals,hold) else hold<-rep("NA",0)
		#if(pass) hold[8] <=0.3 & hold[8] >0) print(c(gA,gB,hold[8]))
		if(pass)l<-l+1
	}
	Pfp[k]<-(l-1)
	k<-k+1
}
print("Complete")
elapsed<-(proc.time()-start)[3]
print(elapsed)


row.names(AllVals)<-c()
AllEpiData<-data.frame(AllVals)

AllEpiData$geneA<-g1
AllEpiData$geneB<-g2
#colnames(AllEpiData)<-c("Bias","Total","Sex","Diet","Bias:gBT","Bias:Sex","Total:Sex","Bias:Diet","Total:Diet","Sex:Diet","Bias:Total:Sex","Bias:Total:Diet","Bias:Sex:Diet","Total:Sex:Diet","Bias:Total:Sex:Diet","adjusted.R.squared","geneA","geneB")

#colnames(AllEpiData)<-c("InterceptP","SexP","DietP","BiasP","TotalP","IntVar","SexVar","DietVar","BiasVar","TotalVar","Intercept","Sex","Diet","Bias","Total","R^2","A","B")


#colnames(AllEpiData)<-c("IntCoe","TotalCoe","Bias:TotalCoe","Total:MHCoe","Total:MLCoe","Total:FHCoe","Total:FLCoe","BT:MHCoe","BT:MLCoe","BT:FHCoe","BT:FLCoe","MeanTErr","VarTErr","gA","gB")

colnames(AllEpiData)<-c("IntCoe","TotalCoe","Bias:TotalCoe","Total:FdCoe","Total:MdCoe","BT:FdCoe","BT:MdCoe","MeanTErr","VarTErr","gA","gB")

###Adjust P-vals


#print("FDR Corrections...")

#Bf<-NULL
#Tf<-NULL
#Sf<-NULL
#Df<-NULL
#BTf<-NULL
#BSf<-NULL
#TSf<-NULL
#BDf<-NULL
#TDf<-NULL
#SDf<-NULL
#BTSf<-NULL
#BTDf<-NULL
#BSDf<-NULL
#TSDf<-NULL
#BTSDf<-NULL

#for(i in 1:length(Pfp)){
#	Bf[Pip[i]:Pfp[i]]<-p.adjust(AllEpiData[Pip[i]:Pfp[i],1],method="fdr")
#	Tf[Pip[i]:Pfp[i]]<-p.adjust(AllEpiData[Pip[i]:Pfp[i],2],method="fdr")
#	Sf[Pip[i]:Pfp[i]]<-p.adjust(AllEpiData[Pip[i]:Pfp[i],3],method="fdr")
#	Df[Pip[i]:Pfp[i]]<-p.adjust(AllEpiData[Pip[i]:Pfp[i],4],method="fdr")
#	BTf[Pip[i]:Pfp[i]]<-p.adjust(AllEpiData[Pip[i]:Pfp[i],5],method="fdr")

#	BSf[Pip[i]:Pfp[i]]<-p.adjust(AllEpiData[Pip[i]:Pfp[i],6],method="fdr")
#	TSf[Pip[i]:Pfp[i]]<-p.adjust(AllEpiData[Pip[i]:Pfp[i],7],method="fdr")
#	BDf[Pip[i]:Pfp[i]]<-p.adjust(AllEpiData[Pip[i]:Pfp[i],8],method="fdr")
#	TDf[Pip[i]:Pfp[i]]<-p.adjust(AllEpiData[Pip[i]:Pfp[i],9],method="fdr")
#	SDf[Pip[i]:Pfp[i]]<-p.adjust(AllEpiData[Pip[i]:Pfp[i],10],method="fdr")
	
#	BTSf[Pip[i]:Pfp[i]]<-p.adjust(AllEpiData[Pip[i]:Pfp[i],11],method="fdr")
#	BTDf[Pip[i]:Pfp[i]]<-p.adjust(AllEpiData[Pip[i]:Pfp[i],12],method="fdr")
#	BSDf[Pip[i]:Pfp[i]]<-p.adjust(AllEpiData[Pip[i]:Pfp[i],13],method="fdr")
#	TSDf[Pip[i]:Pfp[i]]<-p.adjust(AllEpiData[Pip[i]:Pfp[i],14],method="fdr")
#	BTSDf[Pip[i]:Pfp[i]]<-p.adjust(AllEpiData[Pip[i]:Pfp[i],15],method="fdr")
#}
#print("Complete")

#print("Compiling, Filtering, and Exporting")

#AllEpiData$Bias.FDR<-Bf
#AllEpiData$Total.FDR<-Tf
#AllEpiData$Sex.FDR<-Sf
#AllEpiData$Diet.FDR<-Df
#AllEpiData$"Bias:gBT.FDR"<-BTf
#AllEpiData$"Bias:Sex.FDR"<-BSf
#AllEpiData$"Total:Sex.FDR"<-TSf
#AllEpiData$"Bias:Diet.FDR"<-BDf
#AllEpiData$"Total:Diet.FDR"<-TDf
#AllEpiData$"Sex:Diet.FDR"<-SDf
#AllEpiData$"Bias:Total:Sex.FDR"<-BTSf
#AllEpiData$"Bias:Total:Diet.FDR"<-BTDf
#AllEpiData$"Bias:Sex:Diet.FDR"<-BSDf
#AllEpiData$"Total:Sex:Diet.FDR"<-TSDf
#AllEpiData$"Bias:Total:Sex:Diet.FDR"<-BTSDf

#hist(AllEpiData$"Bias",xlim=c(0,1))
#hist(AllEpiData$"Bias.FDR",xlim=c(0,1))


write.table(AllEpiData, file=paste(args[3],args[4],args[5],args[6],args[7],'AllEpiData.tsv',sep="_"), quote=FALSE, sep='\t', row.names=FALSE, col.names=FALSE)

print("Complete")



