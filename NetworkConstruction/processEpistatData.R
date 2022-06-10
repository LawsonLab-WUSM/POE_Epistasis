

rm(list=ls())
set.seed(0)

library(ggplot2)
library(cowplot)
library(MASS)
library(viridis)







plotPermFDRQC<-function(PVALUES,QVALUES,QCUTOFFVALUE,NULL.PVALUES){

	
DF<-data.frame(PVALUES,QVALUES)




QvP<-ggplot(DF,aes(x= PVALUES, y= QVALUES))+
geom_line(size=1)+
xlab("p-values")+
ylab("q-values")+
geom_hline(yintercept=QCUTOFFVALUE,color="red")

DF2<-data.frame(
QvalueCutoff=seq(0,1,0.01),
SigTests=unlist(lapply(seq(0,1,0.01),function(CUTOFF) sum(DF$QVALUES<=CUTOFF,na.rm=TRUE) )),
FalsePositives=unlist(lapply(seq(0,1,0.01),function(CUTOFF) sum(DF$QVALUES<=CUTOFF,na.rm=TRUE) ))*seq(0,1,0.01)
)


# BREAKS<-seq(min(DF2$SigTests),max(DF2$SigTests),round(abs(max(DF2$SigTests)-min(DF2$SigTests))/10,0))

BREAKS<-NULL
INIT<-min(DF2$SigTests,na.rm=TRUE)+0.5
FIN<-max(DF2$SigTests,na.rm=TRUE)


while(INIT<=FIN){
	BREAKS <-c(BREAKS,round(INIT,0))
	INIT<-INIT*2
}


SigvQ<-ggplot(DF2,aes(x=QvalueCutoff, y=SigTests+0.1))+
geom_line(size=1)+
scale_y_continuous(name="Significant tests",trans=scales::pseudo_log_trans(base=2), breaks=c(BREAKS) )+
scale_x_continuous(name="q-value cutoff",breaks=seq(0,1,0.1) )+
geom_vline(xintercept=QCUTOFFVALUE,color="red")






CORRESPONDINGROW<-c(1:length(seq(0,1,0.01)))[seq(0,1,0.01)==QCUTOFFVALUE]

SigvFalse<-ggplot(DF2,aes(x=SigTests,y=FalsePositives))+
geom_line(size=1)+
scale_y_continuous(name="Expected false positives",trans=scales::pseudo_log_trans(base=2), breaks=c(BREAKS) )+
scale_x_continuous(name="Significant tests",trans=scales::pseudo_log_trans(base=2), breaks=c(BREAKS))+
geom_vline(xintercept=DF2[CORRESPONDINGROW,2],color="red")+
theme(axis.text.x = element_text(angle = 45, hjust = 1))





COMB.DF<-data.frame(
pvalue=c(NULL.PVALUES, PVALUES),
Group=c(rep("Null",length(NULL.PVALUES)),rep("Real",length(PVALUES)))
)


# BINWIDTH<-max(DF$PVALUES[DF$QVALUES<=QCUTOFFVALUE],na.rm=TRUE)*10
BINWIDTH<-QCUTOFFVALUE/10

print(BINWIDTH)

COMPARE.DISTS<-ggplot(COMB.DF,aes(pvalue))+
geom_histogram(data=subset(COMB.DF,Group=="Real" & !is.na(pvalue) ),aes(y=..count../sum(..count..)),fill="blue",alpha=0.6,binwidth= BINWIDTH)+
geom_histogram(data=subset(COMB.DF,Group=="Null" & !is.na(pvalue)),aes(y=..count../sum(..count..)),fill="orange",alpha=0.6,binwidth= BINWIDTH)+
scale_x_continuous(expand = c(0, 0))+
scale_y_continuous(expand = c(0, 0))+
ylab("P(pvalue)")


COMBINED<-plot_grid(COMPARE.DISTS, QvP, SigvQ, SigvFalse, ncol=2)

return(COMBINED)

	
}





get_density <- function(x, y, ...) {
  dens <- MASS::kde2d(x, y, ...)
  ix <- findInterval(x, dens$x)
  iy <- findInterval(y, dens$y)
  ii <- cbind(ix, iy)
  return(dens$z[ii])
}




#Import Data
NullepiData<-read.table("Permutation/Perm_merged.tsv",header=FALSE)
colnames(NullepiData)<-c("IntCoe","TotalCoe","Bias_TotalCoe","Total_FdCoe","Total_MdCoe","BT_FdCoe","BT_MdCoe","MeanTErr","VarTErr","Pvalue","gA","gB")

NullepiData<-NullepiData[!is.na(NullepiData$MeanTErr),]

EpiData<-read.table("GLM/merged.tsv",header=FALSE)
colnames(EpiData)<-c("IntCoe","TotalCoe","Bias_TotalCoe","Total_FdCoe","Total_MdCoe","BT_FdCoe","BT_MdCoe","MeanTErr","VarTErr","Pvalue","gA","gB")

EpiData<-EpiData[!is.na(EpiData$MeanTErr),]


apply(EpiData[,2:10],2,function(z) quantile(z[!is.na(z)],prob=c(0.001,0.999)))
Iteration<-1:500

ALPHA<-0.05
Miniteration<-5
IterationWindow<-5
ALPHAerror<-1*10^-4
AccumulatedPVals<-NULL
AlphaQuantileValues<-NULL
AlphaQuantileVars<-NULL
Alpha.sd<-NULL


for (ITER in Iteration){
	
	AccumulatedPVals<-c(AccumulatedPVals ,read.table(print(paste(c("Permutation/Iteration",ITER,"merged.tsv"),collapse="_")),header=FALSE)[,10])
	
	AlphaQuantileValues <-rbind(AlphaQuantileValues,quantile(AccumulatedPVals,type=4,probs= ALPHA,na.rm=TRUE))
	
	if(ITER >= Miniteration){

	Alpha.sd<-sd(AlphaQuantileValues[(length(AlphaQuantileValues)-(IterationWindow-1)):length(AlphaQuantileValues)])
	AlphaQuantileVars<-c(AlphaQuantileVars, Alpha.sd)
	
	}

}



IterationsPlot<-ggplot(data.frame(Iteration=(IterationWindow):(length(AlphaQuantileVars)+(IterationWindow-1)), sd.alpha=AlphaQuantileVars),aes(x=Iteration,y=sd.alpha))+
geom_hline(yintercept=ALPHAerror,color="red",size=0.75)+
scale_y_continuous(labels = function(x) format(x, scientific = TRUE))+
geom_point(size=0.25)+
ylab( expression("p"[alpha]["="][0.05]))

ggsave("IterationsPlot.tiff", IterationsPlot, units="in", width=3, height=2.5)



library(qvalue)


Qval_Perm<-qvalue(NullepiData$Pvalue)

hist(Qval_Perm)
plot(Qval_Perm)


hist(EpiData$Pvalue)

Qval_Real<-qvalue(EpiData$Pvalue)

### Don't hard code this unless the package won't work.
Qval_Real<-qvalue(EpiData$Pvalue,pi0 = 0.715)
#Qval_Real<-qvalue(EpiData$Pvalue,pi0 = 0.5)

EpiData$FDR<-Qval_Real$qvalues


RealDist<-ecdf(EpiData$Pvalue[!is.na(EpiData$Pvalue)])
NullDist<-ecdf(NullepiData$Pvalue[!is.na(NullepiData$Pvalue)])



EpiData$PermutP<-unlist(lapply(EpiData$Pvalue,function(t) (NullDist(t)+1e-16)/(RealDist(t)+1e-16) ))





FDRThreshold<-0.1

PermumationMultipleTestsPlot<-plotPermFDRQC(EpiData$Pvalue, EpiData$PermutP, FDRThreshold, NullepiData$Pvalue)
ggsave("PermumationMultipleTests.tiff", PermumationMultipleTestsPlot)






png(filename="Real_Pvalue_Dist.png", width=500,height=500,type=c("cairo"))
hist(Qval_Real)
dev.off()

png(filename="MultipleTestsPlot.png", width=500,height=500,type=c("cairo"))
plot(Qval_Real)
dev.off()





library(lsr)
NullDist_Pvalue<-ecdf(NullepiData$Pvalue)

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



plot(EpiData$Pvalue,EpiData$FDR,pch=16,cex=0.1,xlim=c(0,1),ylim=c(0,1))



EpiData$density <- get_density(-log10(EpiData$FDR+1e-16), -log10(EpiData$MeanTErr+1e-16), n = 100)


# SignificanceThreshold<-0.005

NetworkVolcanoFDR<-ggplot(EpiData,aes(y=-log10(FDR +1e-16),x= -log10(MeanTErr+1e-16), color = density))+
scale_color_viridis()+
geom_point(size=0.5)+
geom_hline(yintercept= -log10(FDRThreshold),col='grey')+
geom_vline(xintercept= -log10(MTEcutoff+1e-16), col='grey')+
ylab("-Log10( FDR )")+
xlab("-Log10( MTE )")


NetworkVolcanoPermut<-ggplot(subset(EpiData, PermutP!=0),aes(y=-log10(PermutP +1e-16),x= -log10(MeanTErr+1e-16), color = density))+
scale_color_viridis()+
geom_point(size=0.5)+
geom_hline(yintercept= -log10(FDRThreshold),col='grey')+
geom_vline(xintercept= -log10(MTEcutoff+1e-16), col='grey')+
ylab("-Log10( FDR )")+
xlab("-Log10( MTE )")



ggsave("NetworkVolcanoFDR.tiff", NetworkVolcanoFDR)


ggsave("NetworkVolcanoPermut.tiff", NetworkVolcanoPermut)








# filtEpi<-subset(EpiData,Qval<= SignificanceThreshold & MeanTErr <=MTEcutoff)

filtEpi_FDR<-subset(EpiData, FDR <= FDRThreshold & MeanTErr <=MTEcutoff)

#filtEpi_FDR<-subset(EpiData, PermutP <= FDRThreshold & MeanTErr <=MTEcutoff)


# nrow(filtEpi)
nrow(filtEpi_FDR)

# write.table(filtEpi, file='filtEpi.txt', quote=FALSE, sep='\t', row.names=FALSE,col.names=TRUE)
write.table(filtEpi_FDR, file='filtEpi_FDR.txt', quote=FALSE, sep='\t', row.names=FALSE,col.names=TRUE)

write.table(EpiData, file='All_Epi_FDR.txt', quote=FALSE, sep='\t', row.names=FALSE,col.names=TRUE)


# #Percent of A genes passed
# 100*length(unique(filtEpi$gA))/length(unique(EpiData$gA))
# #Percent of B genes Passed
# 100*length(unique(filtEpi$gB))/length(unique(EpiData$gB))
# #Percent of interactions Passed
# 100*nrow(filtEpi)/nrow(EpiData)




#Percent of A genes passed
100*length(unique(filtEpi_FDR$gA))/length(unique(EpiData$gA))
#Percent of B genes Passed
100*length(unique(filtEpi_FDR$gB))/length(unique(EpiData$gB))
#Percent of interactions Passed
100*nrow(filtEpi_FDR)/nrow(EpiData)





BINWIDTH<-FDRThreshold/10

DATDF<-data.frame(
pvalue=c(NullepiData$Pvalue+10^-16,EpiData$Pvalue+10^-16,EpiData$FDR+10^-16),
Group=c(rep("Null",nrow(NullepiData)),rep("Real",nrow(EpiData)),rep("Adjusted",nrow(EpiData)))
)





pval_dist<-ggplot(DATDF,aes(pvalue))+
geom_histogram(data=subset(DATDF,Group=="Real" & !is.na(pvalue) ),aes(y=..count../sum(..count..)),fill="blue",alpha=0.6,binwidth= BINWIDTH)+
geom_histogram(data=subset(DATDF,Group=="Null" & !is.na(pvalue)),aes(y=..count../sum(..count..)),fill="orange",alpha=0.6,binwidth= BINWIDTH)+
ylab("P(pvalue)")


ggsave("ComparingPvalueDist.tiff", pval_dist)




