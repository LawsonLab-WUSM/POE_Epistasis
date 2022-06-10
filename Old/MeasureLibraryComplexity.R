rm(list=ls())
set.seed(0)
library(VGAM)

require("ggplot2")
require("gplots")
require("cowplot")



print("Importing Data...")
Data<-read.table("ProcessedData.txt",header=FALSE)
colnames(Data)<-c("Gene","Cross","Sex","Diet","Number","BC","Lbiase","L","S")
Data$Sbiase<-(1-Data$Lbiase)
Data<-subset(Data,Lbiase != "NA") ## Remove NA samples
Data$Lbiase<-Data$L/(Data$L+Data$S)
Data$Sbiase<-1-Data$Lbiase
print("Complete")
Data$Total<-Data$L+Data$S
Data<-subset(Data, Total>=20) # Require at least 20 reads


## Fitting Function
fitbetabin<-function(N,SampleQuantiles,MaxValue){
	Deviation<-NULL
	Alphas<-NULL
	Beta<-NULL
	Index<-1
	for(i in 1:MaxValue){
		for(z in 1:MaxValue){
			TEMP<-quantile(rbetabinom.ab(N,100,shape1=i,shape2=z)/100,probs=seq(0,1,0.01))
			Deviation[Index]<-as.double(sum(abs((TEMP-SampleQuantiles))))
			Alphas[Index]<-i
			Beta[Index]<-z
			TEMP<-NULL
			Index<-Index+1
		}
	}

	#heatmap.2(matrix(data=Deviation,ncol=MaxValue),trace="none",Colv=FALSE,Rowv=FALSE,key.xlab="Deviation")
	fittedAlpha<-Alphas[match(min(as.double(Deviation)),as.double(Deviation))]
	fittedBeta<-Beta[match(min(as.double(Deviation)),as.double(Deviation))]
	dipersionParameter<-1/(1+fittedAlpha+fittedBeta)
	return(c(fittedAlpha,fittedBeta,dipersionParameter))
}
##

### Get rid of extreme values that prevent fitting
plotData<-Data
Data<-subset(Data, Lbiase<0.95 & Lbiase>0.05) 

LibraryParams<-NULL

for(Bt in as.character(unique(Data$BC))){	
	## QQ fitting
	N<-length(subset(Data, BC == Bt)$Lbiase)
	Stuff<-""
	quanTest<-quantile(subset(Data, BC == Bt)$Lbiase,probs=seq(0,1,0.01))

	Parameters<-fitbetabin(N,quanTest,40)

	bbFit<-quantile(rbetabinom.ab(N,100,shape1=Parameters[1],shape2=Parameters[2])/100,probs=seq(0,1,0.01))

	#plot(bbFit,quanTest,xlim=c(0,1),pch=16,cex=0.5,ylab="Data Quantiles",xlab="Fitted Beta-Binom Quantiles")
	#lines(x=c(0,1),y=c(0,1),col="blue")

	png(filename=paste("LibComp_",Bt,".png",sep=""), width=500,height=500,type=c("cairo"))

	plot(table(round(subset(plotData, BC == Bt)$Lbiase,2)),xlim=c(0,1),main=Bt,xlab="Reference Allele Ratio",ylab="Frequency",col="lightblue")

	K<-1000000

	lines(
	((table(round(rbetabinom.ab(K,100,shape1=Parameters[1],shape2=Parameters[2])/100,2))/K)*N),
	col="orange",
	type='l'
	)

	text(x=0.7,y=round(max(table(round(subset(Data, BC == Bt)$Lbiase,3)))*0.75,0), labels = round(Parameters[3],3))
	dev.off()
	
	LibraryParams <- rbind( data.frame(unique(subset(Data, BC == Bt)[,c(2,3,4,5,6)]),round(Parameters[3],3)), LibraryParams)
	
}


colnames(LibraryParams)[6]<-c("rho")

png(filename=paste("Rho_Acros_Libs.png",sep=""), width=1000,height=600, res=150,type=c("cairo"))

ggplot(data= LibraryParams, aes(x=BC, y=rho)) +
ggtitle("Library complexity - WAT") +
geom_bar(stat="identity", fill=rgb(48/255,108/255,142/255)) +
theme(axis.text.x = element_text(angle=90, vjust=0.2, hjust=1, size=10)) +
ylim(0,round(max(max(LibraryParams$rho),0.05)*1.15,2)) +
xlab("Library") +
ylab(expression(rho)) +
geom_hline(yintercept=0.05,linetype="dashed", color = rgb(172/255,51/255,120/255))

dev.off()





