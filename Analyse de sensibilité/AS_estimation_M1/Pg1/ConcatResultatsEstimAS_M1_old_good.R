setwd("/Users/evergu/Projets_divers/ANR_2016/ABIM/SimulationsAnalyses/AS_estimation_M1")


#Concaténer les résultats
sensires=read.table("res/sensi_0001.tsv",header=T,sep="\t")
for (i in 2:256)
{
    sensires.tmp=read.table(paste('res/sensi_',sprintf("%04d",i),'.tsv',sep=''),header=T,sep="\t")
    sensires=rbind(sensires,sensires.tmp)
}


sensires_2=read.table("res2/sensi_0001.tsv",header=T,sep="\t")
for (i in 2:256)
{
    sensires.tmp=read.table(paste('res2/sensi_',sprintf("%04d",i),'.tsv',sep=''),header=T,sep="\t")
    sensires_2=rbind(sensires_2,sensires.tmp)
}

write.table(sensires_2,paste('sensires_2','.tsv',sep=''),row.names=F,sep="\t")


############################################################
#Deux fichiers pour l'analyse de sensibilité
library(multisensi)
setwd("/Users/evergu/Projets_divers/ANR_2016/ABIM/SimulationsAnalyses/AS_estimation_M1")

sensires2=read.table("sensires_2.tsv",header=TRUE)
sensires2ave=subset(sensires2,irep=="average")
sensires2rep=subset(sensires2,irep!="average")


#myvars <- names(sensires2ave) %in% c("mu_lambda","mu_ch","sd2_lambda","sd2_ch")
#ASdesignInfM1ave=sensires2ave[myvars]

ASdesignM12=expand.grid(m_lambda=c(1,5,10,20), cv_lambda=c(0.25,0.5,0.75,1), m_H=c(0.5,2,5,10), cv_H=c(0.25,0.5,0.75,1))

myvars2 <- names(sensires2ave) %in% c("mu_lambda","mu_ch","sd2_lambda","sd2_ch","sigma","irep")
ASoutputInfM1ave=sensires2ave[!myvars2]
ASoutputInfM1ave$RMSE_mu_lambda=ASoutputInfM1ave$bias_mu_lambda^2+ASoutputInfM1ave$std_mu_lambda^2
ASoutputInfM1ave$RMSE_mu_ch=ASoutputInfM1ave$bias_mu_ch^2+ASoutputInfM1ave$std_mu_ch^2
ASoutputInfM1ave$RMSE_sd2_lambda=ASoutputInfM1ave$bias_sd2_lambda^2+ASoutputInfM1ave$std_sd2_lambda^2
ASoutputInfM1ave$RMSE_sd2_ch=ASoutputInfM1ave$bias_sd2_ch^2+ASoutputInfM1ave$std_sd2_ch^2


res_ASuniv=analysis.anoasg(ASoutputInfM1ave, ASdesignM12, nbcomp = 14, sigma.car = NULL, analysis.args = list(formula = 2, keep.outputs = FALSE))

#"bias_mu_lambda"  "bias_mu_ch"      "bias_sd2_lambda" "bias_sd2_ch"     "bias_sigma"      "std_mu_lambda"   "std_mu_ch"       "std_sd2_lambda"  "std_sd2_ch"      "std_sigma"       "RMSE_mu_lambda"  "RMSE_mu_ch"      "RMSE_sd2_lambda" "RMSE_sd2_ch"
resAS=res_ASuniv$SI
myvars_bias <- names(resAS) %in% c("bias_mu_lambda","bias_mu_ch","bias_sd2_lambda","bias_sd2_ch")
resASbias=resAS[myvars_bias]
myvars_std <- names(resAS) %in% c("std_mu_lambda","std_mu_ch","std_sd2_lambda","std_sd2_ch")
resASstd=resAS[myvars_std]
myvars_RMSE <- names(resAS) %in% c("RMSE_mu_lambda","RMSE_mu_ch","RMSE_sd2_lambda","RMSE_sd2_ch")
resASRMSE=resAS[myvars_RMSE]

names(resASbias)[names(resASbias) == "bias_mu_lambda"] <- "bias_m_lambda"
names(resASbias)[names(resASbias) == "bias_mu_ch"] <- "bias_m_H"
names(resASbias)[names(resASbias) == "bias_sd2_lambda"] <- "bias_s2_lambda"
names(resASbias)[names(resASbias) == "bias_sd2_ch"] <- "bias_s2_H"
round(resASbias,3)

names(resASstd)[names(resASstd) == "std_mu_lambda"] <- "s_m_lambda"
names(resASstd)[names(resASstd) == "std_mu_ch"] <- "s_m_H"
names(resASstd)[names(resASstd) == "std_sd2_lambda"] <- "s_s2_lambda"
names(resASstd)[names(resASstd) == "std_sd2_ch"] <- "s_s2_H"
round(resASstd,3)

names(resASRMSE)[names(resASRMSE) == "RMSE_mu_lambda"] <- "RMSE_m_lambda"
names(resASRMSE)[names(resASRMSE) == "RMSE_mu_ch"] <- "RMSE_m_H"
names(resASRMSE)[names(resASRMSE) == "RMSE_sd2_lambda"] <- "RMSE_s2_lambda"
names(resASRMSE)[names(resASRMSE) == "RMSE_sd2_ch"] <- "RMSE_s2_H"
round(resASRMSE,3)


round(res_ASuniv$SI,3)
colSums(res_ASuniv$SI)
colSums(res_ASuniv$SI[1:4,])
round(res_ASuniv$tSI,3)
round(res_ASuniv$iSI,3)


resASmultiv <- multisensi(model=ASoutputInfM1ave, design=ASdesignM12, reduction=NULL, scale=FALSE)
plot(resASmultiv, normalized = TRUE, color = terrain.colors, gsi.plot = FALSE)

resASmultiv.pca <- multisensi(model=ASoutputInfM1ave, design=ASdesignM12, reduction=basis.ACP, scale=FALSE)
plot(resASmultiv.pca, graph = 1)
summary(resASmultiv.pca, digits = 2)
plot(resASmultiv.pca, graph = 2)
#plot(resASmultiv.pca, graph = 3)
