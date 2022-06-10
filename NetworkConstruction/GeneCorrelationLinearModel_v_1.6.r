rm(list=ls())


print("Importing Data...")

#[ASEdata File] [General Data File] [File Range Start ASE] [File Range Stop ASE] [File Range Start DE] [File Range Stop DE]
args <- commandArgs(trailingOnly = TRUE)


#args <- c("ProcessedData.txt","DE_WAT_Merged_Counts.txt",1,52,1,30)


library(MASS)
library(glmnet)
library(methods)
library(boot)

### POE
Data<-read.table(args[1],header=FALSE)
colnames(Data)<-c("Gene","Cross","Sex","Diet","Number","BC","Lbiase","L","S")
Data<-subset(Data,Lbiase != "NA") ## Remove NA samples
Data$Lbiase<-Data$L/(Data$L+Data$S)
Data$Total<-Data$L+Data$S
Data<-subset(Data, BC !="CCGGACC")



### All
General<-read.table(args[2],header=FALSE)
colnames(General)<-c("Gene","Cross","Sex","Diet","Number","BC","Lbiase","L","S")
General<-subset(General,Lbiase != "NA") ## Remove NA samples
General$Lbiase<-General$L/(General$L+General$S)
General$Total<-General$L+General$S
General<-subset(General, BC !="CCGGACC")


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
print(c(unique(as.character(setA$Gene)),unique(as.character(setB$Gene))))

for(z in 1:permNumber){
        playMod<-NULL
		playMod<-glm(Y~Total+BT+Total:Fd+Total:Md+BT:Fd+BT:Md,dat,family=gaussian)

        Coefs<-rbind(Coefs,playMod$coefficients)
        cv.err<-c(cv.err,cv.glm(dat,playMod,K=10)$delta[2])
}


MODEL<-glm(Y~Total+BT+Total:Fd+Total:Md+BT:Fd+BT:Md,dat,family=gaussian)
MODEL0<-glm(Y~1,dat,family=gaussian)


MODEL.PVALUE<-pchisq(deviance(MODEL0)-deviance(MODEL),df.residual(MODEL0)-df.residual(MODEL),lower.tail=FALSE)


#print(glm(Y~Total+BT+Total:Fd+Total:Md+BT:Fd+BT:Md,dat,family=gaussian))
#print(cv.err)

mCVE<-round(mean(cv.err),2)
vCVE<-round(sd(cv.err),2)
mCoefs<-apply(Coefs,2,mean)



ACVs<-c(mCoefs,mCVE,vCVE,MODEL.PVALUE)

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
		#if(pass) print(hold)
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





colnames(AllEpiData)<-c("IntCoe","TotalCoe","Bias:TotalCoe","Total:FdCoe","Total:MdCoe","BT:FdCoe","BT:MdCoe","MeanTErr","VarTErr","Pvalue","gA","gB")



write.table(AllEpiData, file=paste(args[3],args[4],args[5],args[6],'AllEpiData.tsv',sep="_"), quote=FALSE, sep='\t', row.names=FALSE, col.names=FALSE)

print("Complete")












