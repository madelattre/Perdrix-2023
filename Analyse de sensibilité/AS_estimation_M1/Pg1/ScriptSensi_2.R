#rm(list=ls())
#library(multisensi)
isetparam=as.numeric(commandArgs(trailingOnly=TRUE)[1])

## Tuning de l'algo
niter <- 400
nburnin <- 200
M <- 10
algo.tuning <- list(niter=niter,nburnin=nburnin,M=M)

## Design de l'expérience
#ASdesign_M1=expand.grid(esp_lambda=c(0.1,5,10,100), var_lambda=c(0.25^2,0.5^2,2^2,4^2), esp_ch=c(0.5,2,5,10), var_ch=c(0.1^2,0.25^2,0.5^2,1^2))
ASdesign_M1_2=read.table('ASdesign_M1_2.csv',header=T,sep="\t")

delta <- 1
#D <- c(1.1,2,3,5,10,15,25,50,100,200,400)
D <- c(1.1,2,3,5,7,10,15,25,50,100,200)
N <- 100
## Paramètres du modèle
sigma <- 0.05

#esp_lambda <- 8
#var_lambda <- (0.25*esp_lambda)^2

esp_lambda <-ASdesign_M1_2$esp_lambda[isetparam]
var_lambda <-ASdesign_M1_2$var_lambda[isetparam]
mu_lambda <- log(esp_lambda)-0.5*log(1+var_lambda/((esp_lambda)^2))
sd2_lambda <- log(1+var_lambda/((esp_lambda)^2))

#esp_ch <- 0.5
#var_ch <- (0.25*esp_ch)^2

esp_ch <-ASdesign_M1_2$esp_ch[isetparam]
var_ch <-ASdesign_M1_2$var_ch[isetparam]
mu_ch <- log(esp_ch)-0.5*log(1+var_ch/((esp_ch)^2))
sd2_ch <- log(1+var_ch/((esp_ch)^2))

# true.param contient les moyennes et variances des Gaussiennes
true.param <- list(mu_lambda=mu_lambda,mu_ch=mu_ch,sd2_lambda=sd2_lambda,sd2_ch=sd2_ch,sigma=sigma)

nbrep <- 10

source('Sensi.R')

temp<-Sensi(true.param,N,delta,nbrep,D,algo.tuning)

temp.reformat=cbind(matrix(as.vector(unlist(temp$true.param)),ncol=5,nrow=nbrep,byrow=T),temp$biais,temp$stdFIM)

colnames(temp.reformat)=c(names(temp$true.param),paste('bias',names(temp$true.param),sep='_'),paste('std',names(temp$true.param),sep='_'))

temp.reformat=rbind(temp.reformat,colMeans(temp.reformat))
temp.reformat=as.data.frame(temp.reformat)
temp.reformat$irep=c(paste('rep',1:nbrep,sep='.'),'average')


write.table(temp.reformat,paste('res2/sensi_',sprintf("%04d",isetparam),'.tsv',sep=''),row.names=F,sep="\t")
