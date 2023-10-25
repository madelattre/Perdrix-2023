### Algorithme SAEM

## La variance est donnée par la variance résiduelle seule

SAEM <- function(niter,nburnin,M,data,paraminit,paramfix){
  ## niter     : nombre total d'itérations de l'algorithme
  ## nburnin   : nombre d'itérations pour la phase de burn-in de l'algorithme (nburnin < niter)
  ## paraminit : initialisation des paramètres à estimer
  ## paramfix  : indique si les paramètres sont fixes (1) ou aléatoires (0)
  ## M         : nombre d'itérations du Metropolis-Hastings
  ## data      : dataframe contenant les données, doit contenir une colonne ID identifiant et R réponse fonctionnelle
  
  ## N         : nombre de perdrix/sujets
  
  N <- max(data$ID)
  
  Ntot <- length(data$R)
  
  mul <- rep(NA,niter+1) # lambda
  muH <- rep(NA,niter+1) # ch
  sigma2l <- rep(NA,niter+1) # lambda 
  sigma2H <- rep(NA,niter+1) # ch
  sigmaRes <- rep(NA,niter+1) # écart-type résiduel
  
  # Simulations des paramètres individuels non observés (gaussiens?)
  zl <- rep(NA,N) # lambda 
  zH <- rep(NA,N) # ch
  
  # Statistiques exhaustives pour la matrice d'information de Fisher
  
  diffmat <- matrix(0,N,length(seq.y))
  statexhi <- array(0,dim=c(5,N,niter))
  
  # Suite de pas pour l'approximation stochastique
  
  gamma <- rep(1,niter)
  for (k in (nburnin+1):niter){
    gamma[k] <- 1/(k-nburnin)^(2/3)
  }
  
  # Pour l'estimation des effets fixes
  
  alpha <- 10^(-6/(niter-nburnin))
  
  # Pour améliorer l'exploration de l'espace des paramètres pendant le burnin
  
  kappa <- 1
  
  # Initialisation des paramètres et des paramètres individuels
  
  mul[1] <- paraminit$mul # lambda
  muH[1] <- paraminit$muH # ch
  sigma2l[1] <- paraminit$sigma2l # lambda
  sigma2H[1] <- paraminit$sigma2H # ch
  sigmaRes[1] <- paraminit$sigmaRes # écart-type résiduel
  
  llambda <- rnorm(N,mean=mul[1],sd=sqrt(sigma2l[1]))
  lH      <- rnorm(N,mean=muH[1],sd=sqrt(sigma2H[1]))
  
  s3i <- c()
  
  for (i in 1:N){
    
    zl[i] <- llambda[i]
    zH[i] <- lH[i]
    
    li <- exp(zl[i])
    Hi <- exp(zH[i])
    
    for (m in 2:M){
      
      ### tirage de lambda
      
      phitildel <- rnorm(1,mean=mul[1],sd=sqrt(sigma2l[1]))
      u <- runif(1)
      mi_num <- EspR(exp(phitildel), Hi, seq.y)
      sdi_num <- sigmaRes[1] #sqrt(VarR(exp(phitildel), Hi, seq.y, delta))
      mi_denom <- EspR(li, Hi, seq.y)
      sdi_denom <- sigmaRes[1] #sqrt(VarR(li, Hi, seq.y, delta))
      # taux d'acceptation 
      ind <- which(data$ID==i)
      lograte <- sum(log(dnorm(data$R[ind],mi_num,sdi_num))) - sum(log(dnorm(data$R[ind],mi_denom,sdi_denom)))
      cond <- (log(u)<=min(0,lograte))
      #}
      if (cond){
        zl[i] <- phitildel
        li <- exp(zl[i])
      } 
      
      ### tirage de H
      
      phitildeH <- rnorm(1,mean=muH[1],sd=sqrt(sigma2H[1]))
      v <- runif(1)
      mi_num <- EspR(li, exp(phitildeH), seq.y)
      sdi_num <- sigmaRes[1] #sqrt(VarR(li, exp(phitildeH), seq.y, delta))
      mi_denom <- EspR(li, Hi, seq.y)
      sdi_denom <- sigmaRes[1] #sqrt(VarR(li, Hi, seq.y, delta))
      # taux d'acceptation 
      ind <- which(data$ID==i)
      lograte <- sum(log(dnorm(data$R[ind],mi_num,sdi_num))) - sum(log(dnorm(data$R[ind],mi_denom,sdi_denom)))
      cond <- (log(v)<=min(0,lograte))
      #}
      if (cond){
        zH[i] <- phitildeH
        Hi <- exp(zH[i])
      }
      
    }
    
    ### 
    
    diffmat[i,] <- (data$R[ind]-EspR(li,Hi,seq.y))^2
    # s3i <- c(s3i,sum(data$R[ind]-EspR(li,Hi,seq.y))^2)
  }
  
  # Mise à jour des statistiques S1 et S2
  s1l <- sum(zl)
  s1H <- sum(zH)
  
  s2l <- sum(zl^2)
  s2H <- sum(zH^2)
  
  s3 <- sum(diffmat)
  
  # Statistiques exhaustives pour l'estimation de la FIM
  
  ajustdonnee <- apply(diffmat,1,sum)
  
  STATEXHI<- cbind(zl, zH, zl^2, zH^2, ajustdonnee)
  for (i in 1:N){
    statexhi[,i,1] <- STATEXHI[i,]
  }
   
  
  # Mise à jour du paramètre theta 
  mul[2] <- s1l/N
  muH[2] <- s1H/N
  
  sigma2l[2] <- s2l/N - mul[2]^2
  sigma2H[2] <- s2H/N - muH[2]^2
  
  sigmaRes[2] <- sqrt(s3/Ntot)
  
  for (k in 2:nburnin){
    
    s3i <- c()
    
    for (i in 1:N){
      
      li <- exp(zl[i])
      Hi <- exp(zH[i])
      
      for (m in 2:M){
        
        ### tirage de lambda
        
        phitildel <- rnorm(1,mean=mul[k],sd=sqrt(kappa*sigma2l[k]))
        u <- runif(1)
        mi_num <- EspR(exp(phitildel), Hi, seq.y)
        sdi_num <- sigmaRes[k] #sqrt(VarR(exp(phitildel), Hi, seq.y, delta))
        
        mi_denom <- EspR(li, Hi, seq.y)
        sdi_denom <- sigmaRes[k] #sqrt(VarR(li, Hi, seq.y, delta))
        # taux d'acceptation 
        ind <- which(data$ID==i)
        lograte <- sum(log(dnorm(data$R[ind],mi_num,sdi_num))) - sum(log(dnorm(data$R[ind],mi_denom,sdi_denom))) 
        cond <- (log(u)<=min(0,lograte))
        
        if (cond){
          zl[i] <- phitildel
          li <- exp(zl[i])
        } 
        
        ### tirage de H
        
        phitildeH <- rnorm(1,mean=muH[k],sd=sqrt(kappa*sigma2H[k]))
        v <- runif(1)
        mi_num <- EspR(li, exp(phitildeH), seq.y)
        sdi_num <- sigmaRes[k] #sqrt(VarR(li, exp(phitildeH), seq.y, delta))
        mi_denom <- EspR(li, Hi, seq.y)
        sdi_denom <- sigmaRes[k] #sqrt(VarR(li, Hi, seq.y, delta))
        # taux d'acceptation 
        ind <- which(data$ID==i)
        lograte <- sum(log(dnorm(data$R[ind],mi_num,sdi_num))) - sum(log(dnorm(data$R[ind],mi_denom,sdi_denom)))
        cond <- (log(v)<=min(0,lograte))
        
        if (cond){
          zH[i] <- phitildeH
          Hi <- exp(zH[i])
        }
        
        
      }
      
      diffmat[i,] <- (data$R[ind]-EspR(li,Hi,seq.y))^2
      # s3i <- c(s3i,sum(data$R[ind]-EspR(li,Hi,seq.y))^2)
      
    }
    
    # Mise à jour des statistiques S1, S2 et S3
    s1l <- (1-gamma[k])*s1l + gamma[k]*sum(zl)
    s1H <- (1-gamma[k])*s1H + gamma[k]*sum(zH)
    s2l <- (1-gamma[k])*s2l + gamma[k]*sum(zl^2)
    s2H <- (1-gamma[k])*s2H + gamma[k]*sum(zH^2)
    s3  <- (1-gamma[k])*s3  + gamma[k]*sum(diffmat)
    
    # Statistiques exhaustives pour l'estimation de la FIM
    
    ajustdonnee <- apply(diffmat,1,sum)
    
    STATEXHI<- cbind(zl, zH, zl^2, zH^2, ajustdonnee)
    for (i in 1:N){
      statexhi[,i,k] <- statexhi[,i,k-1]*(1-gamma[k])+gamma[k]*STATEXHI[i,]
    }
    
    
    # Mise à jour du paramètre theta
    mul[k+1] <- s1l/N
    muH[k+1] <- s1H/N
    sigma2l[k+1] <- s2l/N - mul[k+1]^2
    sigma2H[k+1] <- s2H/N - muH[k+1]^2
    sigmaRes[k+1] <- sqrt(s3/Ntot)
  }
  
  for (k in (nburnin+1):niter){
    
    s3i <- c()
    
    for (i in 1:N){
      
      li <- exp(zl[i])
      Hi <- exp(zH[i])
      
      for (m in 2:M){
        
        ### tirage de lambda
        
        phitildel <- rnorm(1,mean=mul[k],sd=sqrt(kappa*sigma2l[k]))
        u <- runif(1)
        mi_num <- EspR(exp(phitildel), Hi, seq.y)
        sdi_num <- sigmaRes[k] #sqrt(VarR(exp(phitildel), Hi, seq.y, delta))
        mi_denom <- EspR(li, Hi, seq.y)
        sdi_denom <- sigmaRes[k] #sqrt(VarR(li, Hi, seq.y, delta))
        # taux d'acceptation 
        ind <- which(data$ID==i)
        lograte <- sum(log(dnorm(data$R[ind],mi_num,sdi_num))) - sum(log(dnorm(data$R[ind],mi_denom,sdi_denom))) 
        cond <- (log(u)<=min(0,lograte))
        
        if (cond){
          zl[i] <- phitildel
          li <- exp(zl[i])
        } 
        
        ### tirage de H
        
        phitildeH <- rnorm(1,mean=muH[k],sd=sqrt(kappa*sigma2H[k]))
        v <- runif(1)
        mi_num <- EspR(li, exp(phitildeH), seq.y)
        sdi_num <- sigmaRes[k] #sqrt(VarR(li, exp(phitildeH), seq.y, delta))
        mi_denom <- EspR(li, Hi, seq.y)
        sdi_denom <- sigmaRes[k] #sqrt(VarR(li, Hi, seq.y, delta))
        # taux d'acceptation 
        ind <- which(data$ID==i)
        lograte <- sum(log(dnorm(data$R[ind],mi_num,sdi_num))) - sum(log(dnorm(data$R[ind],mi_denom,sdi_denom)))
        cond <- (log(v)<=min(0,lograte))
        
        if (cond){
          zH[i] <- phitildeH
          Hi <- exp(zH[i])
        }
        
      }
      
      diffmat[i,] <- (data$R[ind]-EspR(li,Hi,seq.y))^2
      #s3i <- c(s3i,sum(data$R[ind]-EspR(li,Hi,seq.y))^2)
      
    }
    
    # Mise à jour des statistiques S1 et S2
    s1l <- (1-gamma[k])*s1l + gamma[k]*sum(zl)
    s1H <- (1-gamma[k])*s1H + gamma[k]*sum(zH)
    s2l <- (1-gamma[k])*s2l + gamma[k]*sum(zl^2)
    s2H <- (1-gamma[k])*s2H + gamma[k]*sum(zH^2)
    s3  <- (1-gamma[k])*s3  + gamma[k]*sum(diffmat)
    
    # Mise à jour du paramètre theta
    mul[k+1] <- s1l/N
    muH[k+1] <- s1H/N
    sigma2l[k+1] <- (s2l/N - mul[k+1]^2)*(paramfix$l==0) + alpha*sigma2l[k]*(paramfix$l==1)
    sigma2H[k+1] <- (s2H/N - muH[k+1]^2)*(paramfix$H==0) + alpha*sigma2H[k]*(paramfix$H==1)
    sigmaRes[k+1] <- sqrt(s3/Ntot)
    
    # Statistiques exhaustives pour l'estimation de la FIM
    
    ajustdonnee <- apply(diffmat,1,sum)
    
    STATEXHI<- cbind(zl, zH, zl^2, zH^2, ajustdonnee)
    for (i in 1:N){
      statexhi[,i,k] <- statexhi[,i,k-1]*(1-gamma[k])+gamma[k]*STATEXHI[i,]
    }
    
    
  }
  
  
  # Estimation de la matrice d'information de Fisher
  
  deltaindi2 <- matrix(0,5,N)

  deltaindi2[1,]<- (statexhi[1,,niter]-mul[niter+1])/sigma2l[niter+1]
  deltaindi2[2,]<- (statexhi[2,,niter]-muH[niter+1])/sigma2H[niter+1]
  deltaindi2[3,]<- -1/2/sigma2l[niter+1] +1/2/sigma2l[niter+1]^2*(statexhi[3,,niter] - 2*statexhi[1,,niter]*mul[niter+1]+mul[niter+1]^2)
  deltaindi2[4,]<- -1/2/sigma2H[niter+1] +1/2/sigma2H[niter+1]^2*(statexhi[4,,niter] - 2*statexhi[2,,niter]*muH[niter+1]+muH[niter+1]^2)
  deltaindi2[5,]<- -length(seq.y)/2/sigmaRes[niter+1]^2+statexhi[5,,niter]/2/sigmaRes[niter+1]^4

  vnasfinal2 <- matrix(0,5,5)

  for (i in 1:N){
    vnasfinal2 <- vnasfinal2 + deltaindi2[,i]%*%t(deltaindi2[,i])
  }

  vnasfinal2 <- vnasfinal2/N
  
  
  # MCMC pour l'estimation de la log-vraisemblance
  
  nbmcmc <- 10000
  
  p <- matrix(0,N,nbmcmc)
  
  for (i in 1:N){
    for (j in 1:nbmcmc){
      psil <- rnorm(1,mean=mul[niter+1],sd=sqrt(sigma2l[niter+1]))
      psiH <- rnorm(1,mean=muH[niter+1],sd=sqrt(sigma2H[niter+1]))
      mi <- EspR(exp(psil), exp(psiH), seq.y)
      sdi <- sigmaRes[niter+1] #sqrt(VarR(exp(phitildel), Hi, seq.y, delta))
      ind <- which(data$ID==i)
      p[i,j] <- exp(sum(log(dnorm(data$R[ind],mi,sdi))))
    }
  }
  
  # Calcul de la log-vraisemblance
  
  p <- p[,1000:nbmcmc]
  
  llcumsum<-apply(p,1,cumsum)
  for (j in 1:dim(llcumsum)[1]){
    llcumsum[j,] <- log(llcumsum[j,]/j)
  }
  
  ll <- rowSums(llcumsum)
  
  # BIC <- -2*ll[length(ll)] + 4*log(N) + log(N*length(seq.y))
  BIC <- -2*ll[length(ll)] + 5*log(N)
   
  return(list(mul=mul,muH=muH,sigma2H=sigma2H,sigma2l=sigma2l,sigmaRes=sigmaRes,vn=vnasfinal2,ll=ll,BIC=BIC))
  #,BIC=BIC
}

# par(mfrow=c(3,2))
# 
# plot(mul,xlab="Itérations",main="Mu Lambda",type='l')
# abline(a=mu_lambda,b=0,col="red")
# 
# plot(sigmal,xlab="Itérations", main="Sd Lambda",type='l')
# abline(a=sd_lambda,b=0,col="red")
# 
# plot(muH,xlab="Itérations",main="Mu ch",type='l')
# abline(a=mu_ch,b=0,col="red")
# 
# plot(sigmaH,xlab="Itérations", main="Sd ch",type='l')
# abline(a=sd_ch,b=0,col="red")
# 
# plot(sigmaRes,xlab="Itérations", main="Sigma Res",type='l')
# abline(a=sigma,b=0,col="red")

