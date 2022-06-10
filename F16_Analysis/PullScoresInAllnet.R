set.seed(0)
rm(list=ls())

F16SNPsInAllnet<-read.delim("F16_SNPS_within_3Mb_of_WATNet.BED",header=FALSE)
colnames(F16SNPsInAllnet)<-c("Chr","Start","Stop","Type")

ImprintingScores<-read.csv("F16ImprintingScores.csv")

AllnetImprintingScores<-ImprintingScores[apply(ImprintingScores[c(2,5)],1,paste,collapse="_") %in% apply(F16SNPsInAllnet[1:2],1,paste,collapse="_"),]

AdditiveScores<-read.csv("F16AdditiveScores.csv")

AllnetAdditiveScores<-AdditiveScores[apply(AdditiveScores[c(2,5)],1,paste,collapse="_") %in% apply(F16SNPsInAllnet[1:2],1,paste,collapse="_"),]

write.csv(AllnetImprintingScores,"WATNetImprintingScores.csv",quote = FALSE,row.names=FALSE)
write.csv(AllnetAdditiveScores,"WATNetAdditiveScores.csv",quote = FALSE,row.names=FALSE)





