set.seed(0)


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


print("Analyzing Interactions...")

library(lsr)

set.seed(0)

iterations<-1:10

PermQuantiles<-NULL
permMeanQuantiles<-NULL
permMeanDiffQuantiles<-NULL




origData<-Data

for(i in iterations){

##### NULL Model
BcrossP<-NULL
Bcross_SexP<-NULL
Bcross_DietP<-NULL
Bcross_Sex_DietP<-NULL

BcrossF<-NULL
Bcross_SexF<-NULL
Bcross_DietF<-NULL
Bcross_Sex_DietF<-NULL

BcrossE<-NULL
Bcross_SexE<-NULL
Bcross_DietE<-NULL
Bcross_Sex_DietE<-NULL

POE_Cross<-NULL
POE_Male_H<-NULL
POE_Male_L<-NULL
POE_Female_H<-NULL
POE_Female_L<-NULL
POE_Male<-NULL
POE_Female<-NULL
POE_High<-NULL
POE_Low<-NULL

GeneName<-NULL

minCount<-2

genList<-NULL

Data$Cross<-rbinom(nrow(Data),1,0.5) #Make Cross Random
Data$Cross[Data$Cross==0]<-"LXS"
Data$Cross[Data$Cross==1]<-"SXL"

Data$Sex<-rbinom(nrow(Data),1,0.5) #Make Sex Random
Data$Sex[Data$Sex==0]<-"M"
Data$Sex[Data$Sex==1]<-"F"

Data$Diet<-rbinom(nrow(Data),1,0.5) #Make Diet Random
Data$Diet[Data$Diet==0]<-"H"
Data$Diet[Data$Diet==1]<-"L"

l=1
for (gene in unique(Data$Gene)){
	genSubdat<-subset(Data,(Gene==gene))

	At<-FALSE
	St<-FALSE
	Dt<-FALSE
	Ct<-FALSE

	Male_Hf_LXS<-nrow(subset(genSubdat,Cross=="LXS" & Sex=="M" & Diet=="H"))
	Male_Hf_SXL<-nrow(subset(genSubdat,Cross=="SXL" & Sex=="M" & Diet=="H"))
	Male_Lf_LXS<-nrow(subset(genSubdat,Cross=="LXS" & Sex=="M" & Diet=="L"))
	Male_Lf_SXL<-nrow(subset(genSubdat,Cross=="SXL" & Sex=="M" & Diet=="L"))
	Female_Hf_LXS<-nrow(subset(genSubdat,Cross=="LXS" & Sex=="F" & Diet=="H"))
	Female_Hf_SXL<-nrow(subset(genSubdat,Cross=="SXL" & Sex=="F" & Diet=="H"))
	Female_Lf_LXS<-nrow(subset(genSubdat,Cross=="LXS" & Sex=="F" & Diet=="L"))
	Female_Lf_SXL<-nrow(subset(genSubdat,Cross=="SXL" & Sex=="F" & Diet=="L"))
	LXS<-nrow(subset(genSubdat,Cross=="LXS"))
	SXL<-nrow(subset(genSubdat,Cross=="SXL"))

	Male_LXS<-nrow(subset(genSubdat,Cross=="LXS" & Sex=="M"))
	Male_SXL<-nrow(subset(genSubdat,Cross=="SXL" & Sex=="M"))
	Female_LXS<-nrow(subset(genSubdat,Cross=="LXS" & Sex=="F"))
	Female_SXL<-nrow(subset(genSubdat,Cross=="SXL" & Sex=="F"))


	Hf_LXS<-nrow(subset(genSubdat,Cross=="LXS" & Diet=="H"))
	Hf_SXL<-nrow(subset(genSubdat,Cross=="SXL" & Diet=="H"))
	Lf_LXS<-nrow(subset(genSubdat,Cross=="LXS" & Diet=="L"))
	Lf_SXL<-nrow(subset(genSubdat,Cross=="SXL" & Diet=="L"))



	if( Male_Hf_LXS >=minCount & Male_Hf_SXL >=minCount ) POE_Male_H[l]<-mean(subset(genSubdat,Sex=="M" & Diet=="H" & Cross=="SXL")$Lbiase)-mean(subset(genSubdat,Sex=="M" & Diet=="H" & Cross=="LXS")$Lbiase) else POE_Male_H[l]<-"NaN"
	if( Male_Lf_LXS >=minCount & Male_Lf_SXL >=minCount ) POE_Male_L[l]<-mean(subset(genSubdat,Sex=="M" & Diet=="L" & Cross=="SXL")$Lbiase)-mean(subset(genSubdat,Sex=="M" & Diet=="L" & Cross=="LXS")$Lbiase) else POE_Male_L[l]<-"NaN"
	if( Female_Hf_LXS >=minCount & Female_Hf_SXL >=minCount ) POE_Female_H[l]<-mean(subset(genSubdat,Sex=="F" & Diet=="H" & Cross=="SXL")$Lbiase)-mean(subset(genSubdat,Sex=="F" & Diet=="H" & Cross=="LXS")$Lbiase) else POE_Female_H[l]<-"NaN"
	if( Female_Lf_LXS >=minCount & Female_Lf_SXL >=minCount ) POE_Female_L[l]<-mean(subset(genSubdat,Sex=="F" & Diet=="L" & Cross=="SXL")$Lbiase)-mean(subset(genSubdat,Sex=="F" & Diet=="L" & Cross=="LXS")$Lbiase) else POE_Female_L[l]<-"NaN"
	if( Male_LXS >=minCount & Male_SXL >=minCount ) POE_Male[l]<-mean(subset(genSubdat,Sex=="M" & Cross=="SXL")$Lbiase)-mean(subset(genSubdat,Sex=="M" & Cross=="LXS")$Lbiase) else POE_Male[l]<-"NaN"
	if( Female_LXS >=minCount & Female_SXL >=minCount ) POE_Female[l]<-mean(subset(genSubdat,Sex=="F" & Cross=="SXL")$Lbiase)-mean(subset(genSubdat,Sex=="F" & Cross=="LXS")$Lbiase) else POE_Female[l]<-"NaN"
	if( Hf_LXS >=minCount & Hf_SXL >=minCount ) POE_High[l]<-mean(subset(genSubdat, Diet=="H" & Cross=="SXL")$Lbiase)-mean(subset(genSubdat, Diet=="H" & Cross=="LXS")$Lbiase) else POE_High[l]<-"NaN"
	if( Lf_LXS >=minCount & Lf_SXL >=minCount ) POE_Low[l]<-mean(subset(genSubdat, Diet=="L" & Cross=="SXL")$Lbiase)-mean(subset(genSubdat, Diet=="L" & Cross=="LXS")$Lbiase) else POE_Low[l]<-"NaN"
	if( LXS >=minCount & SXL >=minCount ) POE_Cross[l]<-mean(subset(genSubdat, Cross=="SXL")$Lbiase)-mean(subset(genSubdat, Cross=="LXS")$Lbiase) else POE_Cross[l]<-"NaN"



	#All
	if( Male_Hf_LXS >=minCount & Male_Hf_SXL >=minCount & Male_Lf_LXS >=minCount & Male_Lf_SXL >=minCount & Female_Hf_LXS >=minCount & Female_Hf_SXL >=minCount & Female_Lf_LXS >=minCount & Female_Lf_SXL >=minCount) At<-TRUE
	#Sex
	if( Male_Hf_LXS + Male_Lf_LXS >=minCount & Male_Lf_SXL + Male_Hf_SXL >=minCount & Female_Hf_LXS + Female_Lf_LXS >=minCount & Female_Hf_SXL + Female_Lf_SXL >=minCount) St<-TRUE
	#Diet
	if( Male_Hf_LXS + Female_Hf_LXS >=minCount & Male_Hf_SXL + Female_Hf_SXL >=minCount & Male_Lf_LXS + Female_Lf_LXS >=minCount & Male_Lf_SXL + Female_Lf_SXL >=minCount) Dt<-TRUE
	#Cross
	if( LXS >1 & SXL >=minCount) Ct<-TRUE
		

	if( At ) bmod<-aov(Lbiase~Cross+Cross:Sex+Cross:Diet+Cross:Sex:Diet,genSubdat)

	if( !At & St ) bmod<-aov(Lbiase~Cross+Cross:Sex,genSubdat)

	if( !At & !St & Dt ) bmod<-aov(Lbiase~Cross+Cross:Diet,genSubdat)
	
	if( !At & !St & !Dt & Ct ) bmod<-aov(Lbiase~Cross,genSubdat)

	if( !At & !St & !Dt & !Ct ) bmod<-rep("NaN",8)
	
	if( At ) BcrossP[l]<-summary(bmod)[[1]][["Pr(>F)"]][1]
	if( At ) Bcross_SexP[l]<-summary(bmod)[[1]][["Pr(>F)"]][2]
	if( At ) Bcross_DietP[l]<-summary(bmod)[[1]][["Pr(>F)"]][3]
	if( At ) Bcross_Sex_DietP[l]<-summary(bmod)[[1]][["Pr(>F)"]][4]

	if( !At & St ) BcrossP[l]<-summary(bmod)[[1]][["Pr(>F)"]][1]
	if( !At & St ) Bcross_SexP[l]<-summary(bmod)[[1]][["Pr(>F)"]][2]
	if( !At & St ) Bcross_DietP[l]<-"NaN"
	if( !At & St ) Bcross_Sex_DietP[l]<-"NaN"

	if( !At & !St & Dt ) BcrossP[l]<-summary(bmod)[[1]][["Pr(>F)"]][1]
	if( !At & !St & Dt ) Bcross_SexP[l]<-"NaN"
	if( !At & !St & Dt ) Bcross_DietP[l]<-summary(bmod)[[1]][["Pr(>F)"]][2]
	if( !At & !St & Dt ) Bcross_Sex_DietP[l]<-"NaN"

	if( !At & !St & !Dt & Ct ) BcrossP[l]<-summary(bmod)[[1]][["Pr(>F)"]][1]
	if( !At & !St & !Dt & Ct ) Bcross_SexP[l]<-"NaN"
	if( !At & !St & !Dt & Ct ) Bcross_DietP[l]<-"NaN"
	if( !At & !St & !Dt & Ct ) Bcross_Sex_DietP[l]<-"NaN"

	if( !At & !St & !Dt & !Ct ) BcrossP[l]<-"NaN"
	if( !At & !St & !Dt & !Ct ) Bcross_SexP[l]<-"NaN"
	if( !At & !St & !Dt & !Ct ) Bcross_DietP[l]<-"NaN"
	if( !At & !St & !Dt & !Ct ) Bcross_Sex_DietP[l]<-"NaN"

	if( At ) BcrossF[l]<-summary(bmod)[[1]][["F value"]][1]
	if( At ) Bcross_SexF[l]<-summary(bmod)[[1]][["F value"]][2]
	if( At ) Bcross_DietF[l]<-summary(bmod)[[1]][["F value"]][3]
	if( At ) Bcross_Sex_DietF[l]<-summary(bmod)[[1]][["F value"]][4]

	if( !At & St ) BcrossF[l]<-summary(bmod)[[1]][["F value"]][1]
	if( !At & St ) Bcross_SexF[l]<-summary(bmod)[[1]][["F value"]][2]
	if( !At & St ) Bcross_DietF[l]<-"NaN"
	if( !At & St ) Bcross_Sex_DietF[l]<-"NaN"

	if( !At & !St & Dt ) BcrossF[l]<-summary(bmod)[[1]][["F value"]][1]
	if( !At & !St & Dt ) Bcross_SexF[l]<-"NaN"
	if( !At & !St & Dt ) Bcross_DietF[l]<-summary(bmod)[[1]][["F value"]][2]
	if( !At & !St & Dt ) Bcross_Sex_DietF[l]<-"NaN"

	if( !At & !St & !Dt & Ct ) BcrossF[l]<-summary(bmod)[[1]][["F value"]][1]
	if( !At & !St & !Dt & Ct ) Bcross_SexF[l]<-"NaN"
	if( !At & !St & !Dt & Ct ) Bcross_DietF[l]<-"NaN"
	if( !At & !St & !Dt & Ct ) Bcross_Sex_DietF[l]<-"NaN"

	if( !At & !St & !Dt & !Ct ) BcrossF[l]<-"NaN"
	if( !At & !St & !Dt & !Ct ) Bcross_SexF[l]<-"NaN"
	if( !At & !St & !Dt & !Ct ) Bcross_DietF[l]<-"NaN"
	if( !At & !St & !Dt & !Ct ) Bcross_Sex_DietF[l]<-"NaN"


	if( Ct ) bEff<-etaSquared(bmod)[,1]
	
	
	if( At ) BcrossE[l]<-bEff[1]
	if( At ) Bcross_SexE[l]<-bEff[2]
	if( At ) Bcross_DietE[l]<-bEff[3]
	if( At ) Bcross_Sex_DietE[l]<-bEff[4]

	if( !At & St ) BcrossE[l]<-bEff[1]
	if( !At & St ) Bcross_SexE[l]<-bEff[2]
	if( !At & St ) Bcross_DietE[l]<-"NaN"
	if( !At & St ) Bcross_Sex_DietE[l]<-"NaN"


	if( !At & !St & Dt ) BcrossE[l]<-bEff[1]
	if( !At & !St & Dt ) Bcross_SexE[l]<-"NaN"
	if( !At & !St & Dt ) Bcross_DietE[l]<-bEff[2]
	if( !At & !St & Dt ) Bcross_Sex_DietE[l]<-"NaN"


	if( !At & !St & !Dt & Ct ) BcrossE[l]<-bEff[1]
	if( !At & !St & !Dt & Ct ) Bcross_SexE[l]<-"NaN"
	if( !At & !St & !Dt & Ct ) Bcross_DietE[l]<-"NaN"
	if( !At & !St & !Dt & Ct ) Bcross_Sex_DietE[l]<-"NaN"

	if( !At & !St & !Dt & !Ct ) BcrossE[l]<-"NaN"
	if( !At & !St & !Dt & !Ct ) Bcross_SexE[l]<-"NaN"
	if( !At & !St & !Dt & !Ct ) Bcross_DietE[l]<-"NaN"
	if( !At & !St & !Dt & !Ct ) Bcross_Sex_DietE[l]<-"NaN"

	bEff<-NULL

	bmod<-NULL

	l=l+1
}



BcrossP[BcrossP == "Inf"]<-"NaN"
Bcross_SexP[Bcross_SexP == "Inf"]<-"NaN"
Bcross_DietP[Bcross_DietP == "Inf"]<-"NaN"
Bcross_Sex_DietP[Bcross_Sex_DietP == "Inf"]<-"NaN"

BcrossF[BcrossF == "Inf"]<-"NaN"
Bcross_SexF[Bcross_SexF == "Inf"]<-"NaN"
Bcross_DietF[Bcross_DietF == "Inf"]<-"NaN"
Bcross_Sex_DietF[Bcross_Sex_DietF == "Inf"]<-"NaN"

BcrossE[BcrossE == "Inf"]<-"NaN"
Bcross_SexE[Bcross_SexE == "Inf"]<-"NaN"
Bcross_DietE[Bcross_DietE == "Inf"]<-"NaN"
Bcross_Sex_DietE[Bcross_Sex_DietE == "Inf"]<-"NaN"






PermQuantiles<-rbind(PermQuantiles,quantile(as.numeric(Bcross_Sex_DietP),type=4,probs=seq(0,1,0.01),na.rm=TRUE))
permMeanQuantiles<-rbind(permMeanQuantiles,apply(PermQuantiles,2,mean))
permMeanDiffQuantiles<-rbind(permMeanDiffQuantiles,permMeanQuantiles[i-1,]-permMeanQuantiles[i,])

}


plot(y=permMeanDiffQuantiles[,51],iterations[1:length(iterations)-1],type='l',ylab="Quantile Deviation from Mean",xlab="Iteration",ylim=c(-0.01,0.05))



interactions<-data.frame(
	unique(Data$Gene),
	as.numeric(BcrossP),
	as.numeric(Bcross_SexP),
	as.numeric(Bcross_DietP),
	as.numeric(Bcross_Sex_DietP),

	as.numeric(BcrossF),
	as.numeric(Bcross_SexF),
	as.numeric(Bcross_DietF),
	as.numeric(Bcross_Sex_DietF),

	as.numeric(BcrossE),
	as.numeric(Bcross_SexE),
	as.numeric(Bcross_DietE),
	as.numeric(Bcross_Sex_DietE),

	as.numeric(POE_Cross),
	as.numeric(POE_Male_H),
	as.numeric(POE_Male_L),
	as.numeric(POE_Female_H),
	as.numeric(POE_Female_L),
	as.numeric(POE_Male),
	as.numeric(POE_Female),
	as.numeric(POE_High),
	as.numeric(POE_Low)
)

colnames(interactions)<-c(
"Gene",

"B_Cross_P",
"B_Cross_Sex_P",
"B_Cross_Diet_P",
"B_Cross_Sex_Diet_P",

"B_Cross_F",
"B_Cross_Sex_F",
"B_Cross_Diet_F",
"B_Cross_Sex_Diet_F",

"B_cross_E",
"B_Cross_Sex_E",
"B_Cross_Diet_E",
"B_Cross_Sex_Diet_E",

"POE_Cross",
"POE_Male_H",
"POE_Male_L",
"POE_Female_H",
"POE_Female_L",
"POE_Male",
"POE_Female",
"POE_High",
"POE_Low"
)

#########
NullModel<-interactions

Data<-origData




BcrossP<-NULL
Bcross_SexP<-NULL
Bcross_DietP<-NULL
Bcross_Sex_DietP<-NULL

BcrossF<-NULL
Bcross_SexF<-NULL
Bcross_DietF<-NULL
Bcross_Sex_DietF<-NULL

BcrossE<-NULL
Bcross_SexE<-NULL
Bcross_DietE<-NULL
Bcross_Sex_DietE<-NULL

POE_Cross<-NULL
POE_Male_H<-NULL
POE_Male_L<-NULL
POE_Female_H<-NULL
POE_Female_L<-NULL
POE_Male<-NULL
POE_Female<-NULL
POE_High<-NULL
POE_Low<-NULL

l=1
for (gene in unique(Data$Gene)){
	genSubdat<-subset(Data,(Gene==gene))

	At<-FALSE
	St<-FALSE
	Dt<-FALSE
	Ct<-FALSE

	Male_Hf_LXS<-nrow(subset(genSubdat,Cross=="LXS" & Sex=="M" & Diet=="H"))
	Male_Hf_SXL<-nrow(subset(genSubdat,Cross=="SXL" & Sex=="M" & Diet=="H"))
	Male_Lf_LXS<-nrow(subset(genSubdat,Cross=="LXS" & Sex=="M" & Diet=="L"))
	Male_Lf_SXL<-nrow(subset(genSubdat,Cross=="SXL" & Sex=="M" & Diet=="L"))
	Female_Hf_LXS<-nrow(subset(genSubdat,Cross=="LXS" & Sex=="F" & Diet=="H"))
	Female_Hf_SXL<-nrow(subset(genSubdat,Cross=="SXL" & Sex=="F" & Diet=="H"))
	Female_Lf_LXS<-nrow(subset(genSubdat,Cross=="LXS" & Sex=="F" & Diet=="L"))
	Female_Lf_SXL<-nrow(subset(genSubdat,Cross=="SXL" & Sex=="F" & Diet=="L"))
	LXS<-nrow(subset(genSubdat,Cross=="LXS"))
	SXL<-nrow(subset(genSubdat,Cross=="SXL"))

	Male_LXS<-nrow(subset(genSubdat,Cross=="LXS" & Sex=="M"))
	Male_SXL<-nrow(subset(genSubdat,Cross=="SXL" & Sex=="M"))
	Female_LXS<-nrow(subset(genSubdat,Cross=="LXS" & Sex=="F"))
	Female_SXL<-nrow(subset(genSubdat,Cross=="SXL" & Sex=="F"))


	Hf_LXS<-nrow(subset(genSubdat,Cross=="LXS" & Diet=="H"))
	Hf_SXL<-nrow(subset(genSubdat,Cross=="SXL" & Diet=="H"))
	Lf_LXS<-nrow(subset(genSubdat,Cross=="LXS" & Diet=="L"))
	Lf_SXL<-nrow(subset(genSubdat,Cross=="SXL" & Diet=="L"))



	if( Male_Hf_LXS >=minCount & Male_Hf_SXL >=minCount ) POE_Male_H[l]<-mean(subset(genSubdat,Sex=="M" & Diet=="H" & Cross=="SXL")$Lbiase)-mean(subset(genSubdat,Sex=="M" & Diet=="H" & Cross=="LXS")$Lbiase) else POE_Male_H[l]<-"NaN"
	if( Male_Lf_LXS >=minCount & Male_Lf_SXL >=minCount ) POE_Male_L[l]<-mean(subset(genSubdat,Sex=="M" & Diet=="L" & Cross=="SXL")$Lbiase)-mean(subset(genSubdat,Sex=="M" & Diet=="L" & Cross=="LXS")$Lbiase) else POE_Male_L[l]<-"NaN"
	if( Female_Hf_LXS >=minCount & Female_Hf_SXL >=minCount ) POE_Female_H[l]<-mean(subset(genSubdat,Sex=="F" & Diet=="H" & Cross=="SXL")$Lbiase)-mean(subset(genSubdat,Sex=="F" & Diet=="H" & Cross=="LXS")$Lbiase) else POE_Female_H[l]<-"NaN"
	if( Female_Lf_LXS >=minCount & Female_Lf_SXL >=minCount ) POE_Female_L[l]<-mean(subset(genSubdat,Sex=="F" & Diet=="L" & Cross=="SXL")$Lbiase)-mean(subset(genSubdat,Sex=="F" & Diet=="L" & Cross=="LXS")$Lbiase) else POE_Female_L[l]<-"NaN"
	if( Male_LXS >=minCount & Male_SXL >=minCount ) POE_Male[l]<-mean(subset(genSubdat,Sex=="M" & Cross=="SXL")$Lbiase)-mean(subset(genSubdat,Sex=="M" & Cross=="LXS")$Lbiase) else POE_Male[l]<-"NaN"
	if( Female_LXS >=minCount & Female_SXL >=minCount ) POE_Female[l]<-mean(subset(genSubdat,Sex=="F" & Cross=="SXL")$Lbiase)-mean(subset(genSubdat,Sex=="F" & Cross=="LXS")$Lbiase) else POE_Female[l]<-"NaN"
	if( Hf_LXS >=minCount & Hf_SXL >=minCount ) POE_High[l]<-mean(subset(genSubdat, Diet=="H" & Cross=="SXL")$Lbiase)-mean(subset(genSubdat, Diet=="H" & Cross=="LXS")$Lbiase) else POE_High[l]<-"NaN"
	if( Lf_LXS >=minCount & Lf_SXL >=minCount ) POE_Low[l]<-mean(subset(genSubdat, Diet=="L" & Cross=="SXL")$Lbiase)-mean(subset(genSubdat, Diet=="L" & Cross=="LXS")$Lbiase) else POE_Low[l]<-"NaN"
	if( LXS >=minCount & SXL >=minCount ) POE_Cross[l]<-mean(subset(genSubdat, Cross=="SXL")$Lbiase)-mean(subset(genSubdat, Cross=="LXS")$Lbiase) else POE_Cross[l]<-"NaN"



	#All
	if( Male_Hf_LXS >=minCount & Male_Hf_SXL >=minCount & Male_Lf_LXS >=minCount & Male_Lf_SXL >=minCount & Female_Hf_LXS >=minCount & Female_Hf_SXL >=minCount & Female_Lf_LXS >=minCount & Female_Lf_SXL >=minCount) At<-TRUE
	#Sex
	if( Male_Hf_LXS + Male_Lf_LXS >=minCount & Male_Lf_SXL + Male_Hf_SXL >=minCount & Female_Hf_LXS + Female_Lf_LXS >=minCount & Female_Hf_SXL + Female_Lf_SXL >=minCount) St<-TRUE
	#Diet
	if( Male_Hf_LXS + Female_Hf_LXS >=minCount & Male_Hf_SXL + Female_Hf_SXL >=minCount & Male_Lf_LXS + Female_Lf_LXS >=minCount & Male_Lf_SXL + Female_Lf_SXL >=minCount) Dt<-TRUE
	#Cross
	if( LXS >1 & SXL >=minCount) Ct<-TRUE
		

	if( At ) bmod<-aov(Lbiase~Cross+Cross:Sex+Cross:Diet+Cross:Sex:Diet,genSubdat)

	if( !At & St ) bmod<-aov(Lbiase~Cross+Cross:Sex,genSubdat)

	if( !At & !St & Dt ) bmod<-aov(Lbiase~Cross+Cross:Diet,genSubdat)
	
	if( !At & !St & !Dt & Ct ) bmod<-aov(Lbiase~Cross,genSubdat)

	if( !At & !St & !Dt & !Ct ) bmod<-rep("NaN",8)
	
	if( At ) BcrossP[l]<-summary(bmod)[[1]][["Pr(>F)"]][1]
	if( At ) Bcross_SexP[l]<-summary(bmod)[[1]][["Pr(>F)"]][2]
	if( At ) Bcross_DietP[l]<-summary(bmod)[[1]][["Pr(>F)"]][3]
	if( At ) Bcross_Sex_DietP[l]<-summary(bmod)[[1]][["Pr(>F)"]][4]

	if( !At & St ) BcrossP[l]<-summary(bmod)[[1]][["Pr(>F)"]][1]
	if( !At & St ) Bcross_SexP[l]<-summary(bmod)[[1]][["Pr(>F)"]][2]
	if( !At & St ) Bcross_DietP[l]<-"NaN"
	if( !At & St ) Bcross_Sex_DietP[l]<-"NaN"

	if( !At & !St & Dt ) BcrossP[l]<-summary(bmod)[[1]][["Pr(>F)"]][1]
	if( !At & !St & Dt ) Bcross_SexP[l]<-"NaN"
	if( !At & !St & Dt ) Bcross_DietP[l]<-summary(bmod)[[1]][["Pr(>F)"]][2]
	if( !At & !St & Dt ) Bcross_Sex_DietP[l]<-"NaN"

	if( !At & !St & !Dt & Ct ) BcrossP[l]<-summary(bmod)[[1]][["Pr(>F)"]][1]
	if( !At & !St & !Dt & Ct ) Bcross_SexP[l]<-"NaN"
	if( !At & !St & !Dt & Ct ) Bcross_DietP[l]<-"NaN"
	if( !At & !St & !Dt & Ct ) Bcross_Sex_DietP[l]<-"NaN"

	if( !At & !St & !Dt & !Ct ) BcrossP[l]<-"NaN"
	if( !At & !St & !Dt & !Ct ) Bcross_SexP[l]<-"NaN"
	if( !At & !St & !Dt & !Ct ) Bcross_DietP[l]<-"NaN"
	if( !At & !St & !Dt & !Ct ) Bcross_Sex_DietP[l]<-"NaN"

	if( At ) BcrossF[l]<-summary(bmod)[[1]][["F value"]][1]
	if( At ) Bcross_SexF[l]<-summary(bmod)[[1]][["F value"]][2]
	if( At ) Bcross_DietF[l]<-summary(bmod)[[1]][["F value"]][3]
	if( At ) Bcross_Sex_DietF[l]<-summary(bmod)[[1]][["F value"]][4]

	if( !At & St ) BcrossF[l]<-summary(bmod)[[1]][["F value"]][1]
	if( !At & St ) Bcross_SexF[l]<-summary(bmod)[[1]][["F value"]][2]
	if( !At & St ) Bcross_DietF[l]<-"NaN"
	if( !At & St ) Bcross_Sex_DietF[l]<-"NaN"

	if( !At & !St & Dt ) BcrossF[l]<-summary(bmod)[[1]][["F value"]][1]
	if( !At & !St & Dt ) Bcross_SexF[l]<-"NaN"
	if( !At & !St & Dt ) Bcross_DietF[l]<-summary(bmod)[[1]][["F value"]][2]
	if( !At & !St & Dt ) Bcross_Sex_DietF[l]<-"NaN"

	if( !At & !St & !Dt & Ct ) BcrossF[l]<-summary(bmod)[[1]][["F value"]][1]
	if( !At & !St & !Dt & Ct ) Bcross_SexF[l]<-"NaN"
	if( !At & !St & !Dt & Ct ) Bcross_DietF[l]<-"NaN"
	if( !At & !St & !Dt & Ct ) Bcross_Sex_DietF[l]<-"NaN"

	if( !At & !St & !Dt & !Ct ) BcrossF[l]<-"NaN"
	if( !At & !St & !Dt & !Ct ) Bcross_SexF[l]<-"NaN"
	if( !At & !St & !Dt & !Ct ) Bcross_DietF[l]<-"NaN"
	if( !At & !St & !Dt & !Ct ) Bcross_Sex_DietF[l]<-"NaN"

	if( Ct ) bEff<-etaSquared(bmod)[,1]
	
	
	if( At ) BcrossE[l]<-bEff[1]
	if( At ) Bcross_SexE[l]<-bEff[2]
	if( At ) Bcross_DietE[l]<-bEff[3]
	if( At ) Bcross_Sex_DietE[l]<-bEff[4]

	if( !At & St ) BcrossE[l]<-bEff[1]
	if( !At & St ) Bcross_SexE[l]<-bEff[2]
	if( !At & St ) Bcross_DietE[l]<-"NaN"
	if( !At & St ) Bcross_Sex_DietE[l]<-"NaN"


	if( !At & !St & Dt ) BcrossE[l]<-bEff[1]
	if( !At & !St & Dt ) Bcross_SexE[l]<-"NaN"
	if( !At & !St & Dt ) Bcross_DietE[l]<-bEff[2]
	if( !At & !St & Dt ) Bcross_Sex_DietE[l]<-"NaN"


	if( !At & !St & !Dt & Ct ) BcrossE[l]<-bEff[1]
	if( !At & !St & !Dt & Ct ) Bcross_SexE[l]<-"NaN"
	if( !At & !St & !Dt & Ct ) Bcross_DietE[l]<-"NaN"
	if( !At & !St & !Dt & Ct ) Bcross_Sex_DietE[l]<-"NaN"

	if( !At & !St & !Dt & !Ct ) BcrossE[l]<-"NaN"
	if( !At & !St & !Dt & !Ct ) Bcross_SexE[l]<-"NaN"
	if( !At & !St & !Dt & !Ct ) Bcross_DietE[l]<-"NaN"
	if( !At & !St & !Dt & !Ct ) Bcross_Sex_DietE[l]<-"NaN"


	bEff<-NULL

	bmod<-NULL

	l=l+1
}

BcrossP[BcrossP == "Inf"]<-"NaN"
Bcross_SexP[Bcross_SexP == "Inf"]<-"NaN"
Bcross_DietP[Bcross_DietP == "Inf"]<-"NaN"
Bcross_Sex_DietP[Bcross_Sex_DietP == "Inf"]<-"NaN"


BcrossF[BcrossF == "Inf"]<-"NaN"
Bcross_SexF[Bcross_SexF == "Inf"]<-"NaN"
Bcross_DietF[Bcross_DietF == "Inf"]<-"NaN"
Bcross_Sex_DietF[Bcross_Sex_DietF == "Inf"]<-"NaN"


BcrossE[BcrossE == "Inf"]<-"NaN"
Bcross_SexE[Bcross_SexE == "Inf"]<-"NaN"
Bcross_DietE[Bcross_DietE == "Inf"]<-"NaN"
Bcross_Sex_DietE[Bcross_Sex_DietE == "Inf"]<-"NaN"






interactions<-data.frame(
	unique(Data$Gene),
	as.numeric(BcrossP),
	as.numeric(Bcross_SexP),
	as.numeric(Bcross_DietP),
	as.numeric(Bcross_Sex_DietP),

	as.numeric(BcrossF),
	as.numeric(Bcross_SexF),
	as.numeric(Bcross_DietF),
	as.numeric(Bcross_Sex_DietF),

	as.numeric(BcrossE),
	as.numeric(Bcross_SexE),
	as.numeric(Bcross_DietE),
	as.numeric(Bcross_Sex_DietE),

	as.numeric(POE_Cross),
	as.numeric(POE_Male_H),
	as.numeric(POE_Male_L),
	as.numeric(POE_Female_H),
	as.numeric(POE_Female_L),
	as.numeric(POE_Male),
	as.numeric(POE_Female),
	as.numeric(POE_High),
	as.numeric(POE_Low)
)

colnames(interactions)<-c(
"Gene",

"B_Cross_P",
"B_Cross_Sex_P",
"B_Cross_Diet_P",
"B_Cross_Sex_Diet_P",

"B_Cross_F",
"B_Cross_Sex_F",
"B_Cross_Diet_F",
"B_Cross_Sex_Diet_F",

"B_cross_E",
"B_Cross_Sex_E",
"B_Cross_Diet_E",
"B_Cross_Sex_Diet_E",

"POE_Cross",
"POE_Male_H",
"POE_Male_L",
"POE_Female_H",
"POE_Female_L",
"POE_Male",
"POE_Female",
"POE_High",
"POE_Low"
)

png(filename="POE_Score_vs_NULL.png", width=1000,height=500,type=c("cairo"))
par(mfrow=c(1,3),mai = c(0.6, 0.6, 0.1,0.1))

plot(density(NullModel$POE_Cross[!is.na(NullModel$POE_Cross)]),cex.lab=2,lwd=2,col="orange",xlim=c(-1,-0.3),ylim=c(0,2),main="",xlab="")
lines(density(interactions$POE_Cross[!is.na(interactions$POE_Cross)]),lwd=2,col="blue",lty=3)

plot(density(NullModel$POE_Cross[!is.na(NullModel$POE_Cross)]),cex.lab=2,lwd=2,col="orange",xlim=c(-0.3,0.3),ylim=c(0,25),ylab="",main="",xlab="POE Score")
lines(density(interactions$POE_Cross[!is.na(interactions$POE_Cross)]),lwd=2,col="blue",lty=3)

plot(density(NullModel$POE_Cross[!is.na(NullModel$POE_Cross)]),cex.lab=2,lwd=2,col="orange",xlim=c(0.3,1),ylim=c(0,2),ylab="",main="",xlab="")
lines(density(interactions$POE_Cross[!is.na(interactions$POE_Cross)]),lwd=2,col="blue",lty=3)


dev.off()


### Calculating POE score cutoffs
POE_Cross_C<-c(max(abs(quantile(NullModel$POE_Cross[!is.na(NullModel$POE_Cross)],c(0.01,0.99)))))

POE_Male_C<-c(max(abs(quantile(NullModel$POE_Male[!is.na(NullModel$POE_Male)],c(0.01,0.99)))))
POE_Female_C<-c(max(abs(quantile(NullModel$POE_Female[!is.na(NullModel$POE_Female)],c(0.01,0.99)))))

POE_High_C<-c(max(abs(quantile(NullModel$POE_High[!is.na(NullModel$POE_High)],c(0.01,0.99)))))
POE_Low_C<-c(max(abs(quantile(NullModel$POE_Low[!is.na(NullModel$POE_Low)],c(0.01,0.99)))))

POE_Male_H_C<-c(max(abs(quantile(NullModel$POE_Male_H[!is.na(NullModel$POE_Male_H)],c(0.01,0.99)))))
POE_Male_L_C<-c(max(abs(quantile(NullModel$POE_Male_L[!is.na(NullModel$POE_Male_L)],c(0.01,0.99)))))
POE_Female_H_C<-c(max(abs(quantile(NullModel$POE_Female_H[!is.na(NullModel$POE_Female_H)],c(0.01,0.99)))))
POE_Female_L_C<-c(max(abs(quantile(NullModel$POE_Female_L[!is.na(NullModel$POE_Female_L)],c(0.01,0.99)))))


NullDist<-ecdf(NullModel$B_Cross_F[!is.na(NullModel$B_Cross_F)])

interactions$B_Cross_PermutP<-unlist(lapply(interactions$B_Cross_F,function(t) 1-NullDist(t)))
interactions$B_Cross_Sex_PermutP<-unlist(lapply(interactions$B_Cross_Sex_F,function(t) 1-NullDist(t)))
interactions$B_Cross_Diet_PermutP<-unlist(lapply(interactions$B_Cross_Diet_F,function(t) 1-NullDist(t)))
interactions$B_Cross_Sex_Diet_PermutP<-unlist(lapply(interactions$B_Cross_Sex_Diet_F,function(t) 1-NullDist(t)))


interactions$B_Cross_FDR<-p.adjust(interactions$B_Cross_P,method="fdr")
interactions$B_Cross_Sex_FDR<-p.adjust(interactions$B_Cross_Sex_P,method="fdr")
interactions$B_Cross_Diet_FDR<-p.adjust(interactions$B_Cross_Diet_P,method="fdr")
interactions$B_Cross_Sex_Diet_FDR<-p.adjust(interactions$B_Cross_Sex_Diet_P,method="fdr")
write.table(interactions, file='StatsData.tsv', quote=FALSE, sep='\t', row.names=FALSE)


CrossGenes<-subset(interactions, B_Cross_PermutP<=0.05 & abs(POE_Cross)>POE_Cross_C)$Gene
CrossSexGenes<-subset(interactions,B_Cross_Sex_PermutP<=0.05 & ( abs(POE_Male)>POE_Male_C|abs(POE_Female)>POE_Female_C))$Gene
CrossDietGenes<-subset(interactions,B_Cross_Diet_PermutP<=0.05 & ( abs(POE_High)>POE_High_C|abs(POE_Low)>POE_Low_C))$Gene
CrossSexDietGenes<-subset(interactions, B_Cross_Sex_Diet_PermutP<=0.05 & ((abs(POE_Male_H)>POE_Male_H_C | abs(POE_Male_L)>POE_Male_L_C | abs(POE_Female_H)>POE_Female_H_C | abs(POE_Female_L)>POE_Female_L_C)))$Gene


Passed<-union(union(union(CrossGenes,CrossSexGenes),CrossDietGenes),CrossSexDietGenes)



Passed


write.table(Passed, file='Passed.txt', quote=FALSE, sep='\t', row.names=FALSE,col.names=FALSE)

######################
#                    #
#                    #
#      Plotting      #
#                    #
#                    #
######################



makebiasplots<-function(gene){

#Lbias

error.bar <- function(x, y, upper, lower=upper, length=0.1,...){
	if(length(x) != length(y) | length(y) !=length(lower) | length(lower) != length(upper))
	stop("vectors must be same length")
	arrows(x,y+upper, x, y-lower, angle=90, code=3, length=length, ...)
}

crossDf<-cbind(
	if(length(subset(Data,Gene==gene & Cross == "LXS" & Diet == "H" & Sex == "M")$Lbiase)>0) subset(Data,Gene==gene & Cross == "LXS" & Diet == "H" & Sex == "M")$Lbiase else 0,
	if(length(subset(Data,Gene==gene & Cross == "SXL" & Diet == "H" & Sex == "M")$Lbiase)>0) subset(Data,Gene==gene & Cross == "SXL" & Diet == "H" & Sex == "M")$Lbiase else 0,
	if(length(subset(Data,Gene==gene & Cross == "LXS" & Diet == "H" & Sex == "F")$Lbiase)>0) subset(Data,Gene==gene & Cross == "LXS" & Diet == "H" & Sex == "F")$Lbiase else 0,
	if(length(subset(Data,Gene==gene & Cross == "SXL" & Diet == "H" & Sex == "F")$Lbiase)>0) subset(Data,Gene==gene & Cross == "SXL" & Diet == "H" & Sex == "F")$Lbiase else 0,

	if(length(subset(Data,Gene==gene & Cross == "LXS" & Diet == "L" & Sex == "M")$Lbiase)>0) subset(Data,Gene==gene & Cross == "LXS" & Diet == "L" & Sex == "M")$Lbiase else 0,
	if(length(subset(Data,Gene==gene & Cross == "SXL" & Diet == "L" & Sex == "M")$Lbiase)>0) subset(Data,Gene==gene & Cross == "SXL" & Diet == "L" & Sex == "M")$Lbiase else 0,
	if(length(subset(Data,Gene==gene & Cross == "LXS" & Diet == "L" & Sex == "F")$Lbiase)>0) subset(Data,Gene==gene & Cross == "LXS" & Diet == "L" & Sex == "F")$Lbiase else 0,
	if(length(subset(Data,Gene==gene & Cross == "SXL" & Diet == "L" & Sex == "F")$Lbiase)>0) subset(Data,Gene==gene & Cross == "SXL" & Diet == "L" & Sex == "F")$Lbiase else 0
)

crossDf.n<-c(
	length(subset(Data,Gene==gene & Cross == "LXS" & Diet == "H" & Sex == "M")$Lbiase),
	length(subset(Data,Gene==gene & Cross == "SXL" & Diet == "H" & Sex == "M")$Lbiase),
	length(subset(Data,Gene==gene & Cross == "LXS" & Diet == "H" & Sex == "F")$Lbiase),
	length(subset(Data,Gene==gene & Cross == "SXL" & Diet == "H" & Sex == "F")$Lbiase),

	length(subset(Data,Gene==gene & Cross == "LXS" & Diet == "L" & Sex == "M")$Lbiase),
	length(subset(Data,Gene==gene & Cross == "SXL" & Diet == "L" & Sex == "M")$Lbiase),
	length(subset(Data,Gene==gene & Cross == "LXS" & Diet == "L" & Sex == "F")$Lbiase),
	length(subset(Data,Gene==gene & Cross == "SXL" & Diet == "L" & Sex == "F")$Lbiase)
)

colnames(crossDf)<-c("LXS.H.M","SXL.H.M","LXS.H.F","SXL.H.F","LXS.L.M","SXL.L.M","LXS.L.F","SXL.L.F")

crossDf.means<-apply(crossDf,2,mean)
crossDf.sd<-apply(crossDf,2,sd)

png(filename=paste(gene,"_bias_Barplot.png",sep=""), width=500,height=500,type=c("cairo"))
barx <- barplot(crossDf.means, main=gene , ylim=c(0,max(crossDf.means,1/1.3)*1.3), col="blue", axis.lty=1, xlab="Cross", ylab="Lbias")
error.bar(barx,crossDf.means, crossDf.sd/crossDf.n)
text(barx-0.3,crossDf.means+max(crossDf.means)*0.1,crossDf.n,cex=0.7)
dev.off()
}

for(j in union(union(union(CrossGenes,CrossSexGenes),CrossDietGenes),CrossSexDietGenes)){makebiasplots(j)}

