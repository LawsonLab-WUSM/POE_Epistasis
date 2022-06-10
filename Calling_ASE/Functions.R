


runTestsOnGene<-function(DATAFRAME,GENE,makeRandom=FALSE){

	minCount<-2
	genSubdat<-subset(DATAFRAME,(Gene==GENE))


	if(makeRandom==TRUE){
		
		randomIndexes<-sample(1:nrow(genSubdat),nrow(genSubdat))
		
		genSubdat$Cross<-genSubdat$Cross[randomIndexes]
		genSubdat$Sex<-genSubdat$Sex[randomIndexes]
		genSubdat$Diet<-genSubdat$Diet[randomIndexes]
		
	}



	At<-FALSE
	St<-FALSE
	Dt<-FALSE
	Ct<-FALSE

	
	#### Determing counts
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


	#### Calculating group POE scores	
	if( Male_Hf_LXS >=minCount & Male_Hf_SXL >=minCount ) POE_Male_H<-mean(subset(genSubdat,Sex=="M" & Diet=="H" & Cross=="SXL")$Lbiase)-mean(subset(genSubdat,Sex=="M" & Diet=="H" & Cross=="LXS")$Lbiase) else POE_Male_H<-NA
	if( Male_Lf_LXS >=minCount & Male_Lf_SXL >=minCount ) POE_Male_L<-mean(subset(genSubdat,Sex=="M" & Diet=="L" & Cross=="SXL")$Lbiase)-mean(subset(genSubdat,Sex=="M" & Diet=="L" & Cross=="LXS")$Lbiase) else POE_Male_L<-NA
	if( Female_Hf_LXS >=minCount & Female_Hf_SXL >=minCount ) POE_Female_H<-mean(subset(genSubdat,Sex=="F" & Diet=="H" & Cross=="SXL")$Lbiase)-mean(subset(genSubdat,Sex=="F" & Diet=="H" & Cross=="LXS")$Lbiase) else POE_Female_H<-NA
	if( Female_Lf_LXS >=minCount & Female_Lf_SXL >=minCount ) POE_Female_L<-mean(subset(genSubdat,Sex=="F" & Diet=="L" & Cross=="SXL")$Lbiase)-mean(subset(genSubdat,Sex=="F" & Diet=="L" & Cross=="LXS")$Lbiase) else POE_Female_L<-NA
	if( Male_LXS >=minCount & Male_SXL >=minCount ) POE_Male<-mean(subset(genSubdat,Sex=="M" & Cross=="SXL")$Lbiase)-mean(subset(genSubdat,Sex=="M" & Cross=="LXS")$Lbiase) else POE_Male<-NA
	if( Female_LXS >=minCount & Female_SXL >=minCount ) POE_Female<-mean(subset(genSubdat,Sex=="F" & Cross=="SXL")$Lbiase)-mean(subset(genSubdat,Sex=="F" & Cross=="LXS")$Lbiase) else POE_Female<-NA
	if( Hf_LXS >=minCount & Hf_SXL >=minCount ) POE_High<-mean(subset(genSubdat, Diet=="H" & Cross=="SXL")$Lbiase)-mean(subset(genSubdat, Diet=="H" & Cross=="LXS")$Lbiase) else POE_High<-NA
	if( Lf_LXS >=minCount & Lf_SXL >=minCount ) POE_Low<-mean(subset(genSubdat, Diet=="L" & Cross=="SXL")$Lbiase)-mean(subset(genSubdat, Diet=="L" & Cross=="LXS")$Lbiase) else POE_Low<-NA
	if( LXS >=minCount & SXL >=minCount ) POE_Cross<-mean(subset(genSubdat, Cross=="SXL")$Lbiase)-mean(subset(genSubdat, Cross=="LXS")$Lbiase) else POE_Cross<-NA



	



	#### Deciding which tests to run and running them


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

	if( !At & !St & !Dt & !Ct ) bmod<-rep(NA,8)
	
	if( At ) BcrossP<-summary(bmod)[[1]][["Pr(>F)"]][1]
	if( At ) Bcross_SexP<-summary(bmod)[[1]][["Pr(>F)"]][2]
	if( At ) Bcross_DietP<-summary(bmod)[[1]][["Pr(>F)"]][3]
	if( At ) Bcross_Sex_DietP<-summary(bmod)[[1]][["Pr(>F)"]][4]

	if( !At & St ) BcrossP<-summary(bmod)[[1]][["Pr(>F)"]][1]
	if( !At & St ) Bcross_SexP<-summary(bmod)[[1]][["Pr(>F)"]][2]
	if( !At & St ) Bcross_DietP<-NA
	if( !At & St ) Bcross_Sex_DietP<-NA

	if( !At & !St & Dt ) BcrossP<-summary(bmod)[[1]][["Pr(>F)"]][1]
	if( !At & !St & Dt ) Bcross_SexP<-NA
	if( !At & !St & Dt ) Bcross_DietP<-summary(bmod)[[1]][["Pr(>F)"]][2]
	if( !At & !St & Dt ) Bcross_Sex_DietP<-NA

	if( !At & !St & !Dt & Ct ) BcrossP<-summary(bmod)[[1]][["Pr(>F)"]][1]
	if( !At & !St & !Dt & Ct ) Bcross_SexP<-NA
	if( !At & !St & !Dt & Ct ) Bcross_DietP<-NA
	if( !At & !St & !Dt & Ct ) Bcross_Sex_DietP<-NA

	if( !At & !St & !Dt & !Ct ) BcrossP<-NA
	if( !At & !St & !Dt & !Ct ) Bcross_SexP<-NA
	if( !At & !St & !Dt & !Ct ) Bcross_DietP<-NA
	if( !At & !St & !Dt & !Ct ) Bcross_Sex_DietP<-NA

	if( At ) BcrossF<-summary(bmod)[[1]][["F value"]][1]
	if( At ) Bcross_SexF<-summary(bmod)[[1]][["F value"]][2]
	if( At ) Bcross_DietF<-summary(bmod)[[1]][["F value"]][3]
	if( At ) Bcross_Sex_DietF<-summary(bmod)[[1]][["F value"]][4]

	if( !At & St ) BcrossF<-summary(bmod)[[1]][["F value"]][1]
	if( !At & St ) Bcross_SexF<-summary(bmod)[[1]][["F value"]][2]
	if( !At & St ) Bcross_DietF<-NA
	if( !At & St ) Bcross_Sex_DietF<-NA

	if( !At & !St & Dt ) BcrossF<-summary(bmod)[[1]][["F value"]][1]
	if( !At & !St & Dt ) Bcross_SexF<-NA
	if( !At & !St & Dt ) Bcross_DietF<-summary(bmod)[[1]][["F value"]][2]
	if( !At & !St & Dt ) Bcross_Sex_DietF<-NA

	if( !At & !St & !Dt & Ct ) BcrossF<-summary(bmod)[[1]][["F value"]][1]
	if( !At & !St & !Dt & Ct ) Bcross_SexF<-NA
	if( !At & !St & !Dt & Ct ) Bcross_DietF<-NA
	if( !At & !St & !Dt & Ct ) Bcross_Sex_DietF<-NA

	if( !At & !St & !Dt & !Ct ) BcrossF<-NA
	if( !At & !St & !Dt & !Ct ) Bcross_SexF<-NA
	if( !At & !St & !Dt & !Ct ) Bcross_DietF<-NA
	if( !At & !St & !Dt & !Ct ) Bcross_Sex_DietF<-NA


	if( Ct ) bEff<-etaSquared(bmod)[,1]
	
	
	if( At ) BcrossE<-bEff[1]
	if( At ) Bcross_SexE<-bEff[2]
	if( At ) Bcross_DietE<-bEff[3]
	if( At ) Bcross_Sex_DietE<-bEff[4]

	if( !At & St ) BcrossE<-bEff[1]
	if( !At & St ) Bcross_SexE<-bEff[2]
	if( !At & St ) Bcross_DietE<-NA
	if( !At & St ) Bcross_Sex_DietE<-NA


	if( !At & !St & Dt ) BcrossE<-bEff[1]
	if( !At & !St & Dt ) Bcross_SexE<-NA
	if( !At & !St & Dt ) Bcross_DietE<-bEff[2]
	if( !At & !St & Dt ) Bcross_Sex_DietE<-NA


	if( !At & !St & !Dt & Ct ) BcrossE<-bEff[1]
	if( !At & !St & !Dt & Ct ) Bcross_SexE<-NA
	if( !At & !St & !Dt & Ct ) Bcross_DietE<-NA
	if( !At & !St & !Dt & Ct ) Bcross_Sex_DietE<-NA

	if( !At & !St & !Dt & !Ct ) BcrossE<-NA
	if( !At & !St & !Dt & !Ct ) Bcross_SexE<-NA
	if( !At & !St & !Dt & !Ct ) Bcross_DietE<-NA
	if( !At & !St & !Dt & !Ct ) Bcross_Sex_DietE<-NA

	bEff<-NULL

	bmod<-NULL

	c(
	GENE,
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



}
	
	
	

	
	
	
	
	


runPerIndividualPOESCores<-function(DATAFRAME,GENE){

	genSubdat<-subset(DATAFRAME,(Gene==GENE))


	subset(genSubdat,Sex=="M" & Diet=="H" & Cross=="SXL")$BC
	subset(genSubdat,Sex=="M" & Diet=="H" & Cross=="LXS")$BC
	subset(genSubdat,Sex=="M" & Diet=="L" & Cross=="SXL")$BC
	subset(genSubdat,Sex=="M" & Diet=="L" & Cross=="LXS")$BC

	subset(genSubdat,Sex=="F" & Diet=="H" & Cross=="SXL")$BC
	subset(genSubdat,Sex=="F" & Diet=="H" & Cross=="LXS")$BC
	subset(genSubdat,Sex=="F" & Diet=="L" & Cross=="SXL")$BC
	subset(genSubdat,Sex=="F" & Diet=="L" & Cross=="LXS")$BC



		
	M_H_SXL_POESCORE<-subset(genSubdat,Sex=="M" & Diet=="H" & Cross=="SXL")$Lbiase-mean(subset(genSubdat,Sex=="M" & Diet=="H" & Cross=="LXS")$Lbiase)
	
	M_H_LXS_POESCORE<-mean(subset(genSubdat,Sex=="M" & Diet=="H" & Cross=="SXL")$Lbiase)-subset(genSubdat,Sex=="M" & Diet=="H" & Cross=="LXS")$Lbiase
	
	M_L_SXL_POESCORE<-subset(genSubdat,Sex=="M" & Diet=="L" & Cross=="SXL")$Lbiase-mean(subset(genSubdat,Sex=="M" & Diet=="L" & Cross=="LXS")$Lbiase)
	
	M_L_LXS_POESCORE<-mean(subset(genSubdat,Sex=="M" & Diet=="L" & Cross=="SXL")$Lbiase)-subset(genSubdat,Sex=="M" & Diet=="L" & Cross=="LXS")$Lbiase
	
	
	F_H_SXL_POESCORE<-subset(genSubdat,Sex=="F" & Diet=="H" & Cross=="SXL")$Lbiase-mean(subset(genSubdat,Sex=="F" & Diet=="H" & Cross=="LXS")$Lbiase)
	
	F_H_LXS_POESCORE<-mean(subset(genSubdat,Sex=="F" & Diet=="H" & Cross=="SXL")$Lbiase)-subset(genSubdat,Sex=="F" & Diet=="H" & Cross=="LXS")$Lbiase
	
	F_L_SXL_POESCORE<-subset(genSubdat,Sex=="F" & Diet=="L" & Cross=="SXL")$Lbiase-mean(subset(genSubdat,Sex=="F" & Diet=="L" & Cross=="LXS")$Lbiase)
	
	F_L_LXS_POESCORE<-mean(subset(genSubdat,Sex=="F" & Diet=="L" & Cross=="SXL")$Lbiase)-subset(genSubdat,Sex=="F" & Diet=="L" & Cross=="LXS")$Lbiase
	

	data.frame(
	
	
	Gene=c(
	as.character(subset(genSubdat,Sex=="M" & Diet=="H" & Cross=="SXL")$Gene),
	as.character(subset(genSubdat,Sex=="M" & Diet=="H" & Cross=="LXS")$Gene),
	as.character(subset(genSubdat,Sex=="M" & Diet=="L" & Cross=="SXL")$Gene),
	as.character(subset(genSubdat,Sex=="M" & Diet=="L" & Cross=="LXS")$Gene),
	as.character(subset(genSubdat,Sex=="F" & Diet=="H" & Cross=="SXL")$Gene),
	as.character(subset(genSubdat,Sex=="F" & Diet=="H" & Cross=="LXS")$Gene),
	as.character(subset(genSubdat,Sex=="F" & Diet=="L" & Cross=="SXL")$Gene),
	as.character(subset(genSubdat,Sex=="F" & Diet=="L" & Cross=="LXS")$Gene)
	)
	,
POE_Score=c(M_H_SXL_POESCORE,M_H_LXS_POESCORE,M_L_SXL_POESCORE,M_L_LXS_POESCORE,F_H_SXL_POESCORE,F_H_LXS_POESCORE,F_L_SXL_POESCORE,F_L_LXS_POESCORE)
	,
	BC=c(
	as.character(subset(genSubdat,Sex=="M" & Diet=="H" & Cross=="SXL")$BC),
	as.character(subset(genSubdat,Sex=="M" & Diet=="H" & Cross=="LXS")$BC),
	as.character(subset(genSubdat,Sex=="M" & Diet=="L" & Cross=="SXL")$BC),
	as.character(subset(genSubdat,Sex=="M" & Diet=="L" & Cross=="LXS")$BC),
	as.character(subset(genSubdat,Sex=="F" & Diet=="H" & Cross=="SXL")$BC),
	as.character(subset(genSubdat,Sex=="F" & Diet=="H" & Cross=="LXS")$BC),
	as.character(subset(genSubdat,Sex=="F" & Diet=="L" & Cross=="SXL")$BC),
	as.character(subset(genSubdat,Sex=="F" & Diet=="L" & Cross=="LXS")$BC)
	)
	,
	Cross=c(
	as.character(subset(genSubdat,Sex=="M" & Diet=="H" & Cross=="SXL")$Cross),
	as.character(subset(genSubdat,Sex=="M" & Diet=="H" & Cross=="LXS")$Cross),
	as.character(subset(genSubdat,Sex=="M" & Diet=="L" & Cross=="SXL")$Cross),
	as.character(subset(genSubdat,Sex=="M" & Diet=="L" & Cross=="LXS")$Cross),
	as.character(subset(genSubdat,Sex=="F" & Diet=="H" & Cross=="SXL")$Cross),
	as.character(subset(genSubdat,Sex=="F" & Diet=="H" & Cross=="LXS")$Cross),
	as.character(subset(genSubdat,Sex=="F" & Diet=="L" & Cross=="SXL")$Cross),
	as.character(subset(genSubdat,Sex=="F" & Diet=="L" & Cross=="LXS")$Cross)
	)
	,
	Sex=c(
	as.character(subset(genSubdat,Sex=="M" & Diet=="H" & Cross=="SXL")$Sex),
	as.character(subset(genSubdat,Sex=="M" & Diet=="H" & Cross=="LXS")$Sex),
	as.character(subset(genSubdat,Sex=="M" & Diet=="L" & Cross=="SXL")$Sex),
	as.character(subset(genSubdat,Sex=="M" & Diet=="L" & Cross=="LXS")$Sex),
	as.character(subset(genSubdat,Sex=="F" & Diet=="H" & Cross=="SXL")$Sex),
	as.character(subset(genSubdat,Sex=="F" & Diet=="H" & Cross=="LXS")$Sex),
	as.character(subset(genSubdat,Sex=="F" & Diet=="L" & Cross=="SXL")$Sex),
	as.character(subset(genSubdat,Sex=="F" & Diet=="L" & Cross=="LXS")$Sex)
	)
	,
	Diet=c(
	as.character(subset(genSubdat,Sex=="M" & Diet=="H" & Cross=="SXL")$Diet),
	as.character(subset(genSubdat,Sex=="M" & Diet=="H" & Cross=="LXS")$Diet),
	as.character(subset(genSubdat,Sex=="M" & Diet=="L" & Cross=="SXL")$Diet),
	as.character(subset(genSubdat,Sex=="M" & Diet=="L" & Cross=="LXS")$Diet),
	as.character(subset(genSubdat,Sex=="F" & Diet=="H" & Cross=="SXL")$Diet),
	as.character(subset(genSubdat,Sex=="F" & Diet=="H" & Cross=="LXS")$Diet),
	as.character(subset(genSubdat,Sex=="F" & Diet=="L" & Cross=="SXL")$Diet),
	as.character(subset(genSubdat,Sex=="F" & Diet=="L" & Cross=="LXS")$Diet)
	)

	)


}
	
	
	
	
	
	
	
	
summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
    library(plyr)

    # New version of length which can handle NA's: if na.rm==T, don't count them
    length2 <- function (x, na.rm=FALSE) {
        if (na.rm) sum(!is.na(x))
        else       length(x)
    }

    # This does the summary. For each group's data frame, return a vector with
    # N, mean, and sd
    datac <- ddply(data, groupvars, .drop=.drop,
      .fun = function(xx, col) {
        c(N    = length2(xx[[col]], na.rm=na.rm),
          mean = mean   (xx[[col]], na.rm=na.rm),
          sd   = sd     (xx[[col]], na.rm=na.rm)
        )
      },
      measurevar
    )

    # Rename the "mean" column    
    datac <- rename(datac, c("mean" = measurevar))

    datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean

    # Confidence interval multiplier for standard error
    # Calculate t-statistic for confidence interval: 
    # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
    ciMult <- qt(conf.interval/2 + .5, datac$N-1)
    datac$ci <- datac$se * ciMult

    return(datac)
}
	
	
	
	
calcPerGenePOESummary<-function(POE.DF,GENE){
	
	SUMMARY<-summarySE(subset(POE.DF, Gene ==GENE)[,c(2,5,6)],measurevar="POE_Score",groupvars=c("Sex","Diet"))
	
	SUMMARY$Gene<-rep(GENE,nrow(SUMMARY))
	
	return(SUMMARY)
	
}	
	

	
plotPOEScoreBarPlot<-function(POESCORE.DF,GENE){
	
ggplot(summarySE(subset(POESCORE.DF, Gene == GENE)[,c(2,5,6)],measurevar="POE_Score",groupvars=c("Sex","Diet")),aes(x=Sex,y=POE_Score,fill=Diet,color=Sex))+
geom_bar(position=position_dodge(), stat="identity",size=0.5)+
geom_errorbar(aes(ymin= POE_Score-se, ymax= POE_Score +se), width=.2,size=0.5, position=position_dodge(.9))+
ylab("POE Score")+
ylim(-1,1)+
scale_color_manual(values=c("#cb4154","#4169e1"),guide=FALSE)+
scale_fill_manual(values=c("#808080", "#DCDCDC"))+
ggtitle(GENE)

}


	
	
	
	
	
	
	