rm(list=ls())

source('Simus_modM1.R')

## Tuning de l'algo

niter <- 400
nburnin <- 200
M <- 10

## Design de l'expérience

delta <- 1
seq.y <- c(1.1,2,3,5,10,15,25,50,100,200,400)

### Variance résiduelle

sigma <- 0.05

### Simulation des effets aléatoires

esp_lambda <- 8
var_lambda <- (0.25*esp_lambda)^2
mu_lambda <- log(esp_lambda)-0.5*log(1+var_lambda/((esp_lambda)^2))
sd2_lambda <- log(1+var_lambda/((esp_lambda)^2))

esp_ch <- 0.5
var_ch <- (0.25*esp_ch)^2
mu_ch <- log(esp_ch)-0.5*log(1+var_ch/((esp_ch)^2))
sd2_ch <- log(1+var_ch/((esp_ch)^2))


nbsimus <- 100

## Avec bruit de mesure uniquement

source('SAEM_M1_MD_bruit.R')

for (N in c(20,50,100)){
  
  nameFile <- paste('resM1sc2bruit_N',N,'.Rdata',sep="")
  
  estim <- list()
  datasim <- list()
  paraminit <- list()
  
  for (k in 1:nbsimus){
    
    ### Simulation des effets aléatoires
    
    lambda <- exp(rnorm(N,mean=mu_lambda,sd=sqrt(sd2_lambda)))
    ch <- exp(rnorm(N,mean=mu_ch,sd=sqrt(sd2_ch)))
    
    ### Simulation des données de réponse fonctionnelle
    
    R <- rep(0,N*length(seq.y))
    
    for (i in 1:N){
      for (j in 1:length(seq.y)){
        esp <- EspR(lambda[i], ch[i], seq.y[j])
        R[(i-1)*length(seq.y)+j] <- rnorm(1,mean=esp,sd=sigma)
      }
    }
    
    perdrix <- data.frame(ID = rep(seq(1,N,1),each=length(seq.y)), R = R, y = rep(seq.y,times=N))

    paraminit[[k]] <- list(mul=runif(1,mu_lambda-0.5,mu_lambda+0.5),
                      muH=runif(1,mu_ch-0.2,mu_ch+0.2),
                      sigma2l=runif(1,0.5*sd2_lambda,2*sd2_lambda),
                      sigma2H=runif(1,0.5*sd2_ch,2*sd2_ch),
                      sigmaRes=runif(1,1.5*sigma,3*sigma))
    
    paramfix <- list(l=0,H=0)
    
    datasim[[k]] <- perdrix
    
    estim[[k]] <- SAEM(niter,nburnin,M,perdrix,paraminit[[k]],paramfix)
    
    # par(mfrow=c(3,2))
    # 
    # plot(estim[[k]]$mul,xlab="Itérations",main="Mu Lambda",type='l')
    # abline(a=mu_lambda,b=0,col="red")
    # 
    # plot(estim[[k]]$sigma2l,xlab="Itérations", main="Sd2 Lambda",type='l')
    # abline(a=sd2_lambda,b=0,col="red")
    # 
    # plot(estim[[k]]$muH,xlab="Itérations",main="Mu ch",type='l')
    # abline(a=mu_ch,b=0,col="red")
    # 
    # plot(estim[[k]]$sigma2H,xlab="Itérations", main="Sd2 ch",type='l')
    # abline(a=sd2_ch,b=0,col="red")
    # 
    # plot(estim[[k]]$sigmaRes,xlab="Itérations", main="Sigma Res",type='l')
    # abline(a=sigma,b=0,col="red")
    
  }
  
  res <- list(estim=estim,init=paraminit, datasim=datasim)
  
  save(res,file=nameFile)
}



## Avec bruit du modèle et bruit de mesure

rm(list=ls())

source('Simus_modM1.R')

## Tuning de l'algo

niter <- 400
nburnin <- 200
M <- 10

## Design de l'expérience

delta <- 1
seq.y <- c(1.1,2,3,5,10,15,25,50,100,200,400)

### Variance résiduelle

sigma <- 0.05

### Simulation des effets aléatoires

esp_lambda <- 8
var_lambda <- (0.25*esp_lambda)^2
mu_lambda <- log(esp_lambda)-0.5*log(1+var_lambda/((esp_lambda)^2))
sd2_lambda <- log(1+var_lambda/((esp_lambda)^2))

esp_ch <- 0.5
var_ch <- (0.25*esp_ch)^2
mu_ch <- log(esp_ch)-0.5*log(1+var_ch/((esp_ch)^2))
sd2_ch <- log(1+var_ch/((esp_ch)^2))


nbsimus <- 100

source('SAEM_M1_MD_bruit_meca.R')

for (N in c(20,50,100)){
  
  nameFile <- paste('resM1sc2bruitmeca_N',N,'.Rdata',sep="")
  
  estim <- list()
  datasim <- list()
  paraminit <- list()
  
  for (k in 1:nbsimus){
    
    ### Simulation des effets aléatoires
    
    lambda <- exp(rnorm(N,mean=mu_lambda,sd=sqrt(sd2_lambda)))
    ch <- exp(rnorm(N,mean=mu_ch,sd=sqrt(sd2_ch)))
    
    ### Simulation des données de réponse fonctionnelle
    
    R <- rep(0,N*length(seq.y))
    
    for (i in 1:N){
      for (j in 1:length(seq.y)){
        esp <- EspR(lambda[i], ch[i], seq.y[j])
        var <- VarR(lambda[i], ch[i], seq.y[j], delta)
        R[(i-1)*length(seq.y)+j] <- rnorm(1,mean=esp,sd=sqrt(var+sigma^2)) 
      }
    }
    
    
    perdrix <- data.frame(ID = rep(seq(1,N,1),each=length(seq.y)), R = R, D = rep(seq.y,times=N))
    
    datasim[[k]] <- perdrix
    
    paraminit[[k]] <- list(mul=runif(1,mu_lambda-0.5,mu_lambda+0.5),
                      muH=runif(1,mu_ch-0.2,mu_ch+0.2),
                      sigma2l=runif(1,0.5*sd2_lambda,2*sd2_lambda),
                      sigma2H=runif(1,0.5*sd2_ch,2*sd2_ch),
                      sigmaRes=runif(1,1.5*sigma,3*sigma))
    
    paramfix <- list(l=0,H=0)
    
    estim[[k]] <- SAEM(niter,nburnin,M,perdrix,paraminit[[k]],paramfix,delta)
    
    # par(mfrow=c(3,2))
    # 
    # plot(estim[[k]]$mul,xlab="Itérations",main="Mu Lambda",type='l')
    # abline(a=mu_lambda,b=0,col="red")
    # 
    # plot(estim[[k]]$sigma2l,xlab="Itérations", main="Sd2 Lambda",type='l')
    # abline(a=sd2_lambda,b=0,col="red")
    # 
    # plot(estim[[k]]$muH,xlab="Itérations",main="Mu ch",type='l')
    # abline(a=mu_ch,b=0,col="red")
    # 
    # plot(estim[[k]]$sigma2H,xlab="Itérations", main="Sd2 ch",type='l')
    # abline(a=sd2_ch,b=0,col="red")
    # 
    # plot(estim[[k]]$sigmaRes,xlab="Itérations", main="Sigma Res",type='l')
    # abline(a=sigma,b=0,col="red")
    
  }
  
  res <- list(estim=estim,init=paraminit, datasim=datasim)
  
  save(res,file=nameFile)

}


## Avec bruit du modèle uniquement

rm(list=ls())

source('Simus_modM1.R')

## Tuning de l'algo

niter <- 400
nburnin <- 200
M <- 10

## Design de l'expérience

delta <- 1
seq.y <- c(1.1,2,3,5,10,15,25,50,100,200,400)

### Variance résiduelle

sigma <- 0.05

### Simulation des effets aléatoires

esp_lambda <- 8
var_lambda <- (0.25*esp_lambda)^2
mu_lambda <- log(esp_lambda)-0.5*log(1+var_lambda/((esp_lambda)^2))
sd2_lambda <- log(1+var_lambda/((esp_lambda)^2))

esp_ch <- 0.5
var_ch <- (0.25*esp_ch)^2
mu_ch <- log(esp_ch)-0.5*log(1+var_ch/((esp_ch)^2))
sd2_ch <- log(1+var_ch/((esp_ch)^2))


nbsimus <- 100

source('SAEM_M1_MD_meca2.R')

for (N in c(20,50,100)){
  
  nameFile <- paste('resM1sc2meca_N',N,'.Rdata',sep="")
  
  estim <- list()
  datasim <- list()
  paraminit <- list()
  
  for (k in 1:nbsimus){
    
    ### Simulation des effets aléatoires
    
    lambda <- exp(rnorm(N,mean=mu_lambda,sd=sqrt(sd2_lambda)))
    ch <- exp(rnorm(N,mean=mu_ch,sd=sqrt(sd2_ch)))
    
    ### Simulation des données de réponse fonctionnelle
    
    R <- rep(0,N*length(seq.y))
    
    for (i in 1:N){
      for (j in 1:length(seq.y)){
        esp <- EspR(lambda[i], ch[i], seq.y[j])
        var <- VarR(lambda[i], ch[i], seq.y[j], delta)
        R[(i-1)*length(seq.y)+j] <- rnorm(1,mean=esp,sd=sqrt(var)) 
      }
    }
    
    
    perdrix <- data.frame(ID = rep(seq(1,N,1),each=length(seq.y)), R = R, D = rep(seq.y,times=N))
    
    datasim[[k]] <- perdrix
    
    paraminit[[k]] <- list(mul=runif(1,mu_lambda-0.5,mu_lambda+0.5),
                           muH=runif(1,mu_ch-0.2,mu_ch+0.2),
                           sigma2l=runif(1,0.5*sd2_lambda,2*sd2_lambda),
                           sigma2H=runif(1,0.5*sd2_ch,2*sd2_ch),
                           sigmaRes=0.1)
    
    paramfix <- list(l=0,H=0)
    
    estim[[k]] <- SAEM(niter,nburnin,M,perdrix,paraminit[[k]],paramfix,delta)
    
    # par(mfrow=c(2,2))
    # 
    # plot(estim[[k]]$mul,xlab="Itérations",main="Mu Lambda",type='l')
    # abline(a=mu_lambda,b=0,col="red")
    # 
    # plot(estim[[k]]$sigma2l,xlab="Itérations", main="Sd2 Lambda",type='l')
    # abline(a=sd2_lambda,b=0,col="red")
    # 
    # plot(estim[[k]]$muH,xlab="Itérations",main="Mu ch",type='l')
    # abline(a=mu_ch,b=0,col="red")
    # 
    # plot(estim[[k]]$sigma2H,xlab="Itérations", main="Sd2 ch",type='l')
    # abline(a=sd2_ch,b=0,col="red")
    
  }
  
  res <- list(estim=estim,init=paraminit, datasim=datasim)
  
  save(res,file=nameFile)
  
}

###################
#### Résultats ####
###################


library(ggplot2)
library(gridExtra)
library(cowplot)


# 1- modèle avec bruit de mesure uniquement

### Variance résiduelle

sigma <- 0.05

### Simulation des effets aléatoires

esp_lambda <- 8
var_lambda <- (0.25*esp_lambda)^2
mu_lambda <- log(esp_lambda)-0.5*log(1+var_lambda/((esp_lambda)^2))
sd2_lambda <- log(1+var_lambda/((esp_lambda)^2))

esp_ch <- 0.5
var_ch <- (0.25*esp_ch)^2
mu_ch <- log(esp_ch)-0.5*log(1+var_ch/((esp_ch)^2))
sd2_ch <- log(1+var_ch/((esp_ch)^2))


nbsimus <- 100

load('resM1sc2meca_N20.Rdata')

estim <- res$estim

biais <- c()
param <- c()

for (k in 1:nbsimus){
  biais <- c(biais,estim[[k]]$mul[length(estim[[k]]$mul)]-mu_lambda)
  param <- c(param,'muL')
  biais <- c(biais,estim[[k]]$sigma2l[length(estim[[k]]$sigma2l)]-sd2_lambda)
  param <- c(param,'sd2L')
  biais <- c(biais,estim[[k]]$muH[length(estim[[k]]$muH)]-mu_ch)
  param <- c(param,'muH')
  biais <- c(biais,estim[[k]]$sigma2H[length(estim[[k]]$sigma2H)]-sd2_ch)
  param <- c(param,'sd2H')
  # biais <- c(biais,estim[[k]]$sigmaRes[length(estim[[k]]$sigmaRes)]-sigma)
  # param <- c(param,'sigma')
}


res20 <- data.frame(biais=biais,param=param)

p20 <- ggplot(res20, aes(x=param,y=biais, fill=param)) + geom_boxplot() +
  theme(plot.title=element_text(face="bold"),axis.text.x=element_text(size=14),legend.position="none") + xlab("") + ylab("mu") + 
  ggtitle("N=20") + ylim(-0.15,0.15)

load('resM1sc2meca_N50.Rdata')

estim <- res$estim

biais <- c()
param <- c()

for (k in 1:nbsimus){
  biais <- c(biais,estim[[k]]$mul[length(estim[[k]]$mul)]-mu_lambda)
  param <- c(param,'muL')
  biais <- c(biais,estim[[k]]$sigma2l[length(estim[[k]]$sigma2l)]-sd2_lambda)
  param <- c(param,'sd2L')
  biais <- c(biais,estim[[k]]$muH[length(estim[[k]]$muH)]-mu_ch)
  param <- c(param,'muH')
  biais <- c(biais,estim[[k]]$sigma2H[length(estim[[k]]$sigma2H)]-sd2_ch)
  param <- c(param,'sd2H')
  # biais <- c(biais,estim[[k]]$sigmaRes[length(estim[[k]]$sigmaRes)]-sigma)
  # param <- c(param,'sigma')
}


res50 <- data.frame(biais=biais,param=param)

p50 <- ggplot(res50, aes(x=param,y=biais, fill=param)) + geom_boxplot() +
  theme(plot.title=element_text(face="bold"),axis.text.x=element_text(size=14),legend.position="none") + xlab("") + ylab("mu") + 
  ggtitle("N=50") + ylim(-0.15,0.15)

load('resM1sc2meca_N100.Rdata')

estim <- res$estim

biais <- c()
param <- c()

for (k in 1:nbsimus){
  biais <- c(biais,estim[[k]]$mul[length(estim[[k]]$mul)]-mu_lambda)
  param <- c(param,'muL')
  biais <- c(biais,estim[[k]]$sigma2l[length(estim[[k]]$sigma2l)]-sd2_lambda)
  param <- c(param,'sd2L')
  biais <- c(biais,estim[[k]]$muH[length(estim[[k]]$muH)]-mu_ch)
  param <- c(param,'muH')
  biais <- c(biais,estim[[k]]$sigma2H[length(estim[[k]]$sigma2H)]-sd2_ch)
  param <- c(param,'sd2H')
  # biais <- c(biais,estim[[k]]$sigmaRes[length(estim[[k]]$sigmaRes)]-sigma)
  # param <- c(param,'sigma')
}


res100 <- data.frame(biais=biais,param=param)

p100 <- ggplot(res100, aes(x=param,y=biais, fill=param)) + geom_boxplot() +
  theme(plot.title=element_text(face="bold"),axis.text.x=element_text(size=14),legend.position="none") + xlab("") + ylab("mu") + 
  ggtitle("N=100") + ylim(-0.15,0.15) 

plot_grid(p20, p50, p100, ncol = 3, nrow = 1)
