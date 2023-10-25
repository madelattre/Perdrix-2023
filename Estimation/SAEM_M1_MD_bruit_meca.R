### Algorithme SAEM

## La variance est donnée par la variance résiduelle seule

## Remarque : bien gérer le seq.y, comme sur jeu de données réelles et mettre delta en argument


SAEM <- function(niter,nburnin,M,data,paraminit,paramfix,delta){
  ## niter     : nombre total d'itérations de l'algorithme
  ## nburnin   : nombre d'itérations pour la phase de burn-in de l'algorithme (nburnin < niter)
  ## paraminit : initialisation des paramètres à estimer
  ## paramfix  : indique si les paramètres sont fixes (1) ou aléatoires (0)
  ## M         : nombre d'itérations du Metropolis-Hastings
  ## data      : dataframe contenant les données, doit contenir une colonne ID identifiant et R réponse fonctionnelle
  
  ## N         : nombre de perdrix/sujets
  
  N <- max(data$ID)
  
  Ntot <- length(data$R)
  
  seq.y <- unique(data$y)
  
  mul <- rep(NA,niter+1) # lambda
  muH <- rep(NA,niter+1) # ch
  sigma2l <- rep(NA,niter+1) # lambda 
  sigma2H <- rep(NA,niter+1) # ch
  sigmaRes <- rep(NA,niter+1) # écart-type résiduel
  
  # Simulations des paramètres individuels non observés (gaussiens?)
  zl <- rep(NA,N) # lambda 
  zH <- rep(NA,N) # ch
  
  # Suite de pas pour l'approximation stochastique
  
  gamma <- rep(1,niter)
  for (k in (nburnin+1):niter){
    gamma[k] <- 1/(k-nburnin)^(2/3)
  }
  
  # Pour l'estimation des effets fixes
  
  alpha <- 10^(-6/(niter-nburnin))
  
  # Pour améliorer l'exploration de l'espace des paramètres pendant le burnin
  
  kappa <- 1
  
  # Pour Fisher
  
  deltaindi <- array(0,c(5,N,niter+1)) # Dernier indice? 
  tempderiveeas <- matrix(0,5,N)
  
  # Initialisation des paramètres et des paramètres individuels
  
  mul[1] <- paraminit$mul # lambda
  muH[1] <- paraminit$muH # ch
  sigma2l[1] <- paraminit$sigma2l # lambda
  sigma2H[1] <- paraminit$sigma2H # ch
  sigmaRes[1] <- paraminit$sigmaRes # écart-type résiduel
  
  llambda <- rnorm(N,mean=mul[1],sd=sqrt(sigma2l[1])) 
  lH      <- rnorm(N,mean=muH[1],sd=sqrt(sigma2H[1]))
  
  yij <- c()
  moyij <- c()
  varij <- c()
  
  for (i in 1:N){
    
    zl[i] <- llambda[i]
    zH[i] <- lH[i]
    
    li <- exp(zl[i])
    Hi <- exp(zH[i])
    
    for (m in 2:M){
      
      ### tirage de lambda
      
      phitildel <- rnorm(1,mean=mul[1],sd=sqrt(kappa*sigma2l[1]))
      u <- runif(1)
      mi_num <- EspR(exp(phitildel), Hi, seq.y)
      sdi_num <- sqrt(VarR(exp(phitildel), Hi, seq.y, delta)+sigmaRes[1]^2)
      #print(sdi_num)
      mi_denom <- EspR(li, Hi, seq.y)
      sdi_denom <- sqrt(VarR(li, Hi, seq.y, delta)+sigmaRes[1]^2)
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
      
      phitildeH <- rnorm(1,mean=muH[1],sd=sqrt(kappa*sigma2H[1]))
      v <- runif(1)
      mi_num <- EspR(li, exp(phitildeH), seq.y)
      sdi_num <- sqrt(VarR(li, exp(phitildeH), seq.y, delta)+sigmaRes[1]^2)
      mi_denom <- EspR(li, Hi, seq.y)
      sdi_denom <- sqrt(VarR(li, Hi, seq.y, delta)+sigmaRes[1]^2)
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
    
    ## Vecteur des moyennes et vecteur des variances pour toutes les observations
    yij <- c(yij,data$R[ind])
    moyij <- c(moyij, EspR(li, Hi, seq.y))
    varij <- c(varij, VarR(li, Hi, seq.y, delta))

  }
  
  
  # Mise à jour des statistiques S1 et S2
  s1l <- sum(zl)
  s1H <- sum(zH)
  
  s2l <- sum(zl^2)
  s2H <- sum(zH^2)
  
  #resa <- optim(par=sigmaRes[1]^2,f=floglik,y=yij,moyij=moyij,varij=varij)
  resa <- optimize(interval=c(0.0001,5),f=floglik,y=yij,moyij=moyij,varij=varij)
  
  # Mise à jour du paramètre theta 
  mul[2] <- s1l/N
  muH[2] <- s1H/N
  
  sigma2l[2] <- s2l/N - mul[2]^2
  sigma2H[2] <- s2H/N - muH[2]^2
  
  sigmaRes[2] <- sqrt(resa$minimum)
  
  # Pour l'estimation de la matrice d'information de Fisher
  diffmat <- matrix(-1/2*((sigmaRes[2]^2+varij)) + 1/2*(yij-moyij)^2/(varij+sigmaRes[2]^2)^2,nrow=N,byrow=TRUE)
  
  tempderiveeas[1,]<- (zl-mul[2])/sigma2l[2]
  tempderiveeas[2,]<- (zH-muH[2])/sigma2H[2]
  tempderiveeas[3,]<- -1/2/sigma2l[2] +1/2/sigma2l[2]^2*(zl-mul[2])^2
  tempderiveeas[4,]<- -1/2/sigma2H[2] +1/2/sigma2H[2]^2*(zH-muH[2])^2
  tempderiveeas[5,]<- rowSums(diffmat)
  
  deltaindi[,,2] <- deltaindi[,,1]*(1-gamma[2])+gamma[2]*tempderiveeas
  
  for (k in 2:nburnin){
    
    moyij <- c()
    varij <- c()
    
    for (i in 1:N){
      
      li <- exp(zl[i])
      Hi <- exp(zH[i])
      
      for (m in 2:M){
        
        ### tirage de lambda
        
        phitildel <- rnorm(1,mean=mul[k],sd=sqrt(kappa*sigma2l[k]))
        u <- runif(1)
        mi_num <- EspR(exp(phitildel), Hi, seq.y)
        sdi_num <- sqrt(VarR(exp(phitildel), Hi, seq.y, delta)+sigmaRes[k]^2)
        
        mi_denom <- EspR(li, Hi, seq.y)
        sdi_denom <- sqrt(VarR(li, Hi, seq.y, delta)+sigmaRes[k]^2)
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
        sdi_num <- sqrt(VarR(li, exp(phitildeH), seq.y, delta)+sigmaRes[k]^2)
        mi_denom <- EspR(li, Hi, seq.y)
        sdi_denom <- sqrt(VarR(li, Hi, seq.y, delta)+sigmaRes[k]^2)
        # taux d'acceptation 
        ind <- which(data$ID==i)
        lograte <- sum(log(dnorm(data$R[ind],mi_num,sdi_num))) - sum(log(dnorm(data$R[ind],mi_denom,sdi_denom)))
        cond <- (log(v)<=min(0,lograte))
        
        if (cond){
          zH[i] <- phitildeH
          Hi <- exp(zH[i])
        }
        
        
      }
      
      ## Vecteur des moyennes et vecteur des variances pour toutes les observations
      moyij <- c(moyij, EspR(li, Hi, seq.y))
      varij <- c(varij, VarR(li, Hi, seq.y, delta))      
    }
    
    # Mise à jour des statistiques S1, S2
    s1l <- (1-gamma[k])*s1l + gamma[k]*sum(zl)
    s1H <- (1-gamma[k])*s1H + gamma[k]*sum(zH)
    s2l <- (1-gamma[k])*s2l + gamma[k]*sum(zl^2)
    s2H <- (1-gamma[k])*s2H + gamma[k]*sum(zH^2)
    resa <- optimize(interval=c(0.001,5),f=floglik,y=yij,moyij=moyij,varij=varij)
    
    # Mise à jour du paramètre theta
    mul[k+1] <- s1l/N
    muH[k+1] <- s1H/N
    sigma2l[k+1] <- s2l/N - mul[k+1]^2
    sigma2H[k+1] <- s2H/N - muH[k+1]^2
    sigmaRes[k+1] <- sqrt(resa$minimum)
    
    # Pour l'estimation de la matrice d'information de Fisher
    
    diffmat <- matrix(-1/2*((sigmaRes[k+1]^2+varij)) + 1/2*(yij-moyij)^2/(varij+sigmaRes[k+1]^2)^2,nrow=N,byrow=TRUE)
    
    tempderiveeas[1,]<- (zl-mul[k+1])/sigma2l[k+1]
    tempderiveeas[2,]<- (zH-muH[k+1])/sigma2H[k+1]
    tempderiveeas[3,]<- -1/2/sigma2l[k+1] +1/2/sigma2l[k+1]^2*(zl-mul[k+1])^2
    tempderiveeas[4,]<- -1/2/sigma2H[k+1] +1/2/sigma2H[k+1]^2*(zH-muH[k+1])^2
    tempderiveeas[5,]<- rowSums(diffmat)
    
    deltaindi[,,k+1] <- deltaindi[,,k]*(1-gamma[k])+gamma[k]*tempderiveeas
    
  }
  
  for (k in (nburnin+1):niter){
    
    moyij <- c()
    varij <- c()
    
    for (i in 1:N){
      
      li <- exp(zl[i])
      Hi <- exp(zH[i])
      
      for (m in 2:M){
        
        ### tirage de lambda
        
        phitildel <- rnorm(1,mean=mul[k],sd=sqrt(kappa*sigma2l[k]))
        u <- runif(1)
        mi_num <- EspR(exp(phitildel), Hi, seq.y)
        sdi_num <- sqrt(VarR(exp(phitildel), Hi, seq.y, delta)+sigmaRes[k]^2)
        mi_denom <- EspR(li, Hi, seq.y)
        sdi_denom <- sqrt(VarR(li, Hi, seq.y, delta)+sigmaRes[k]^2)
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
        sdi_num <- sqrt(VarR(li, exp(phitildeH), seq.y, delta)+sigmaRes[k]^2)
        mi_denom <- EspR(li, Hi, seq.y)
        sdi_denom <- sqrt(VarR(li, Hi, seq.y, delta)+sigmaRes[k]^2)
        # taux d'acceptation 
        ind <- which(data$ID==i)
        lograte <- sum(log(dnorm(data$R[ind],mi_num,sdi_num))) - sum(log(dnorm(data$R[ind],mi_denom,sdi_denom)))
        cond <- (log(v)<=min(0,lograte))
        
        if (cond){
          zH[i] <- phitildeH
          Hi <- exp(zH[i])
        }
        
      }
      
      ## Vecteur des moyennes et vecteur des variances pour toutes les observations
      moyij <- c(moyij, EspR(li, Hi, seq.y))
      varij <- c(varij, VarR(li, Hi, seq.y, delta))      
      
      
    }
    
    # Mise à jour des statistiques S1 et S2
    s1l <- (1-gamma[k])*s1l + gamma[k]*sum(zl)
    s1H <- (1-gamma[k])*s1H + gamma[k]*sum(zH)
    s2l <- (1-gamma[k])*s2l + gamma[k]*sum(zl^2)
    s2H <- (1-gamma[k])*s2H + gamma[k]*sum(zH^2)
    resa <- optimize(interval=c(0.0001,5),f=floglik,y=yij,moyij=moyij,varij=varij)
    
    # Mise à jour du paramètre theta
    mul[k+1] <- s1l/N
    muH[k+1] <- s1H/N
    # sigmal[k+1] <- sqrt(s2l/N - mul[k+1]^2)*(paramfix$l==0) + sqrt(alpha)*sigmal[k]*(paramfix$l==1)
    # sigmaH[k+1] <- sqrt(s2H/N - muH[k+1]^2)*(paramfix$H==0) + sqrt(alpha)*sigmaH[k]*(paramfix$H==1)
    sigma2l[k+1] <- (s2l/N - mul[k+1]^2)*(paramfix$l==0) + alpha*sigma2l[k]*(paramfix$l==1)
    sigma2H[k+1] <- (s2H/N - muH[k+1]^2)*(paramfix$H==0) + alpha*sigma2H[k]*(paramfix$H==1)
    sigmaRes[k+1] <- (1-gamma[k])*sigmaRes[k] + gamma[k]*sqrt(resa$minimum)
    
    # Pour l'estimation de la matrice d'information de Fisher
    
    diffmat <- matrix(-1/2*((sigmaRes[k+1]^2+varij)) + 1/2*(yij-moyij)^2/(varij+sigmaRes[k+1]^2)^2,nrow=N,byrow=TRUE)
    
    tempderiveeas[1,]<- (zl-mul[k+1])/sigma2l[k+1]
    tempderiveeas[2,]<- (zH-muH[k+1])/sigma2H[k+1]
    tempderiveeas[3,]<- -1/2/sigma2l[k+1] +1/2/sigma2l[k+1]^2*(zl-mul[k+1])^2
    tempderiveeas[4,]<- -1/2/sigma2H[k+1] +1/2/sigma2H[k+1]^2*(zH-muH[k+1])^2
    tempderiveeas[5,]<- rowSums(diffmat)
    
    deltaindi[,,k+1] <- deltaindi[,,k]*(1-gamma[k])+gamma[k]*tempderiveeas
    
  }
  
  # Estimation de la matrice d'information de Fisher
  
  FIM <- matrix(0,5,5)

  FIM <- deltaindi[,,k+1]%*%t(deltaindi[,,k+1])
 
  # MCMC pour l'estimation de la log-vraisemblance
  
  nbmcmc <- 10000
  
  p <- matrix(0,N,nbmcmc)
  
  for (i in 1:N){
    for (j in 1:nbmcmc){
      psil <- rnorm(1,mean=mul[niter+1],sd=sqrt(sigma2l[niter+1]))
      psiH <- rnorm(1,mean=muH[niter+1],sd=sqrt(sigma2H[niter+1]))
      mi <- EspR(exp(psil), exp(psiH), seq.y)
      sdi <- sqrt(VarR(exp(psil), exp(psiH), seq.y, delta)+sigmaRes[niter+1]^2) 
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
  
  #BIC <- -2*ll[length(ll)] + 4*log(N) + log(N*length(seq.y))
  BIC <- -2*ll[length(ll)] + 5*log(N) 
  
  return(list(mul=mul,muH=muH,sigma2H=sigma2H,sigma2l=sigma2l,sigmaRes=sigmaRes,FIM=FIM,ll=ll,BIC=BIC))
}


floglik <- function(a2,y,moyij,varij){
  value <- 1/2*sum(log(a2+varij)) + 1/2*sum((y-moyij)^2/(varij+a2))
  return(value)
}