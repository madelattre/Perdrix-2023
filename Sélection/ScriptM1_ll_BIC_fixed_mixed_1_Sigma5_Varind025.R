## On simule les données sous le modèle où le bruit n'est donné que par la résiduelle

rm(list=ls())

source('Simus_modM1.R')

library(bbmle)

## Tuning de l'algo

niter <- 400
nburnin <- 200
M <- 10

## Design de l'expérience

delta <- 10
seq.y <- c(1.1,2,3,5,10,15,25,50,100,200,400)

### Variance résiduelle

sigma.seq <- c(0.05)

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

N <- 100

for (sigma in sigma.seq){
  
  namefile <- paste('ResBIC_fixed_mixed_M1_1_sigma',sigma*100,'.Rdata',sep="")
  load(namefile)
  estim1 <- res$estim1 
  estim2 <- res$estim2
  estim3 <- res$estim3
  datasim <- res$datasim
  paraminit <- res$paraminit
  ll <- matrix(NA,nbsimus,3)
  BIC <- res$BIC#matrix(NA,nbsimus,6)
  choice <- res$choice#rep(NA,nbsimus)
  coefs1 <- res$coefs1#matrix(NA,nbsimus,3)
  coefs2 <- res$coefs2#matrix(NA,nbsimus,3)
  coefs3 <- res$coefs3#matrix(NA,nbsimus,2)
  
  
  print('Sigma value')
  print(sigma)

  for (k in 98:nbsimus){
    
    print(k)
    
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
                           sigmaRes=runif(1,3*sigma,5*sigma))
    
    paramfix <- list(l=0,H=0)
    
    datasim[[k]] <- perdrix
    
    source('SAEM_M1_MD_bruit.R')
    estim1[[k]] <- SAEM(niter,nburnin,M,perdrix,paraminit[[k]],paramfix)
    source('SAEM_M1_MD_bruit_meca.R')
    estim2[[k]] <- SAEM(niter,nburnin,M,perdrix,paraminit[[k]],paramfix,delta)
    source('SAEM_M1_MD_meca2.R')
    estim3[[k]] <- SAEM(niter,nburnin,M,perdrix,paraminit[[k]],paramfix,delta)
    
    # estimation pour M1.1
    
    minusloglM1_1 <- function(lambda,h,sigma){
      ## - 2 fois log-vraisemblance M1 avec variance résiduelle seule
      yobs <- datasim[[k]]$R # les observations
      esp <- EspR(lambda, h, datasim[[k]]$y)
      ll <- sum(log(sigma^2)+(yobs-esp)^2/sigma^2)
      return(ll)
    }
    
    estim <- mle2(minusloglM1_1, 
                  start = list(lambda = exp(paraminit[[k]]$mul),h = exp(paraminit[[k]]$muH),sigma = paraminit[[k]]$sigmaRes)
    )
    
    coefs1[k,] <- estim@coef
    ll[k,1] <- estim@min
    BIC[k,1] <- ll[k,1] + 3*log(length(datasim[[k]]$R))
    
    # estimation pour M1.2
    
    minusloglM1_2 <- function(lambda,h,sigma){
      ## - 2 fois log-vraisemblance M1 avec variance résiduelle seule
      yobs <- datasim[[k]]$R # les observations
      esp <- EspR(lambda, h, datasim[[k]]$y)
      varr <- VarR(lambda,h, datasim[[k]]$y,delta)
      ll <- sum(log(varr+sigma^2)+(yobs-esp)^2/(varr+sigma^2))
      return(ll)
    }
    
    estim <- mle2(minusloglM1_2, 
                  start = list(lambda = exp(paraminit[[k]]$mul),h = exp(paraminit[[k]]$muH),sigma = paraminit[[k]]$sigmaRes)
    )
    
    coefs2[k,] <- estim@coef
    ll[k,2] <- estim@min
    BIC[k,2] <- ll[k,2] + 3*log(length(datasim[[k]]$R))
    
    
    # estimation pour M1.3
    
    minusloglM1_3 <- function(lambda,h){
      ## - 2 fois log-vraisemblance M1 avec variance résiduelle seule
      yobs <- datasim[[k]]$R # les observations
      esp <- EspR(lambda, h, datasim[[k]]$y)
      varr <- VarR(lambda,h, datasim[[k]]$y,delta)
      ll <- sum(log(varr)+(yobs-esp)^2/(varr))
      return(ll)
    }
    
    estim <-mle2(minusloglM1_3,
                 start = list(lambda = exp(paraminit[[k]]$mul),h = exp(paraminit[[k]]$muH)),
                 method="Nelder-Mead"
    )
    
    
    coefs3[k,] <- estim@coef
    ll[k,3] <- estim@min
    BIC[k,3] <- ll[k,3] + 2*log(length(datasim[[k]]$R))
    
    BIC[k,4] <- estim1[[k]]$BIC
    BIC[k,5] <- estim2[[k]]$BIC
    BIC[k,6] <- estim3[[k]]$BIC
    
    choice[k]<-which.min(BIC[k,])
    
    res <- list(estim1=estim1,estim2=estim2,estim3=estim3,choice=choice,sigma=sigma,coefs1=coefs1,coefs2=coefs2,coefs3=coefs3,BIC=BIC,var.ind=0.25)

    namefile <- paste('ResBIC_fixed_mixed_M1_1_sigma',sigma*100,'_2.Rdata',sep="")

    save(res,file=namefile)
  }
  
  
}
