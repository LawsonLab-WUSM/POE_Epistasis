rm(list=ls())
set.seed(0)
#Import Data
NullepiData<-read.table("Perm_merged.tsv",header=FALSE)
colnames(NullepiData)<-c("IntCoe","TotalCoe","Bias_TotalCoe","Total_FdCoe","Total_MdCoe","BT_FdCoe","BT_MdCoe","MeanTErr","VarTErr","gA","gB")

NullepiData<-NullepiData[!is.na(NullepiData$MeanTErr),]

EpiData<-read.table("merged.tsv",header=FALSE)
colnames(EpiData)<-c("IntCoe","TotalCoe","Bias_TotalCoe","Total_FdCoe","Total_MdCoe","BT_FdCoe","BT_MdCoe","MeanTErr","VarTErr","gA","gB")

EpiData<-EpiData[!is.na(EpiData$MeanTErr),]


apply(EpiData[,2:9],2,function(z) quantile(z[!is.na(z)],prob=c(0.001,0.999)))
Iteration<-1:30

PermQuantiles<-NULL
permMeanQuantiles<-NULL
permMeanDiffQuantiles<-NULL


for (ITER in Iteration){
	PermQuantiles<-rbind(PermQuantiles,quantile(read.table(print(paste(c("Iteration",ITER,"merged.tsv"),collapse="_")),header=FALSE)[,8],type=4,probs=seq(0,1,0.01),na.rm=TRUE))

	permMeanQuantiles<-rbind(permMeanQuantiles,apply(PermQuantiles,2,mean))
	permMeanDiffQuantiles<-rbind(permMeanDiffQuantiles,permMeanQuantiles[ITER-1,]-permMeanQuantiles[ITER,])
}


tiff(filename="WAT_Permutation_mean_deviation.tiff",width = 600, height = 600, res=100)
plot(y=permMeanDiffQuantiles[,51], Iteration[1:length(Iteration)-1],type='l',ylab="Quantile Deviation from Mean",xlab="Iteration")
dev.off()


library(lsr)
NullDist_MTE<-ecdf(NullepiData$MeanTErr)

MTEcutoff<-quantile(NullepiData$MeanTErr,prob=0.01,rm.na=TRUE)


png(filename="NullCumDistPlot.png", width=500,height=500,type=c("cairo"))

plot(
quantile(NullepiData$MeanTErr,prob=seq(0.0001,1,0.001)),
log="x",
y=seq(0.0001,1,0.001),
xlab="MTE",
ylab="P(X<=MTE)",
type='l',
main="Null distribution of MTE"
)

lines(x=c(MTEcutoff,MTEcutoff),c(-1,2),col="red")

dev.off()


plot(y=table(NullepiData$MeanTErr)/sum(table(NullepiData$MeanTErr)),x=as.numeric(rownames(table(NullepiData$MeanTErr))),log="x",type='l',ylim=c(0,0.02))
lines(y=table(EpiData$MeanTErr)/sum(table(EpiData$MeanTErr)),x=as.numeric(rownames(table(EpiData$MeanTErr))),col="blue")



png(filename="CompareCumDistPlot.png", width=500,height=500,type=c("cairo"))

plot(
quantile(NullepiData$MeanTErr,prob=seq(0.0001,1,0.001)),
log="x",
y=seq(0.0001,1,0.001),
xlab="MTE",
ylab="P(X<=MTE)",
type='l',
main="Real versus Null distributions of MTE"
)

lines(
quantile(EpiData$MeanTErr,prob=seq(0.0001,1,0.001)),
y=seq(0.0001,1,0.001),
col="blue"
)

lines(x=c(MTEcutoff,MTEcutoff),c(-1,2),col="red")

dev.off()


EpiData$Pval<-NullDist_MTE(EpiData$MeanTErr)

filtEpi<-subset(EpiData,Pval<=0.01)


write.table(filtEpi, file='filtEpi.txt', quote=FALSE, sep='\t', row.names=FALSE,col.names=TRUE)


#Percent of A genes passed
100*length(unique(filtEpi$gA))/length(unique(EpiData$gA))
#Percent of B genes Passed
100*length(unique(filtEpi$gB))/length(unique(EpiData$gB))
#Percent of interactions Passed
100*nrow(filtEpi)/nrow(EpiData)

nrow(filtEpi)


