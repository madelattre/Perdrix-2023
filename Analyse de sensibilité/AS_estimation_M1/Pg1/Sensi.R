Sensi <- function(true.param,N,delta,nbrep=1,D,algo.tuning){
  #true.param  : liste des paramètres du modèle
  #N           : nb de perdrix
  #delta       : pas de temps lié au modèle mécaniste
  #nbrep       : nb de répétitions par point de grille
  #D           : vecteur des densités de graines
  #algo.tuning : liste des paramètres de tuning de l'algorithme SAEM (nb d'itérations, ...)
  
  biais <- matrix(0,nbrep,5)
  std <- matrix(0,nbrep,5)
  
  source('Simus_modM1.R')
  
  # 1- Simulation d'un jeu de données
  
  mu_lambda <- true.param$mu_lambda
  mu_ch <- true.param$mu_ch
  sd2_lambda <- true.param$sd2_lambda
  sd2_ch <- true.param$sd2_ch
  sigma <- true.param$sigma
  
  for (rep in 1:nbrep){
    
    ## Simulation des effets aléatoires
    
    lambda <- exp(rnorm(N,mean=mu_lambda,sd=sqrt(sd2_lambda)))
    ch <- exp(rnorm(N,mean=mu_ch,sd=sqrt(sd2_ch)))
    
    
    ## Simulation des données de réponse fonctionnelle
    
    R <- rep(0,N*length(D))
    
    for (i in 1:N){
      for (j in 1:length(D)){
        esp <- EspR(lambda[i], ch[i], D[j])
        var <- VarR(lambda[i], ch[i], D[j], delta)
        R[(i-1)*length(D)+j] <- rnorm(1,mean=esp,sd=sqrt(var+sigma^2)) 
        
      }
    }
    
    ## Création d'un data.frame
    
    perdrix <- data.frame(ID = rep(seq(1,N,1),each=length(D)), R = R, y = rep(D,times=N))
    
    # 2- Estimation
    
    source('SAEM_M1_MD_bruit_meca.R')
    
    ## Choix (aléatoire) de valeurs initiales pour l'algorithme
    
    if (mu_lambda>0) {mul_inf=0.9*mu_lambda; mul_sup=1.1*mu_lambda} else {mul_inf=1.1*mu_lambda; mul_sup=0.9*mu_lambda}
    if (mu_ch>0) {much_inf=0.9*mu_ch; much_sup=1.1*mu_ch} else {much_inf=1.1*mu_ch;
        much_sup=0.9*mu_ch}

    paraminit <- list(mul=runif(1,mul_inf,mul_sup), #runif(1,mu_lambda-0.5,mu_lambda+0.5),
                        muH=runif(1,much_inf,much_sup), #runif(1,mu_ch-0.2,mu_ch+0.2),
                      sigma2l=runif(1,0.5*sd2_lambda,2*sd2_lambda),
                      sigma2H=runif(1,0.5*sd2_ch,2*sd2_ch),
                      sigmaRes=0.5)#runif(1,1.5*sigma,3*sigma))
                      
    ## On travaille dans un modèle où les deux paramètres sont aléatoires
    
    paramfix <- list(l=0,H=0)
    
    ## Algo SAEM
    
    niter <- algo.tuning$niter
    nburnin <- algo.tuning$nburnin
    M <- algo.tuning$M
    
    estim <- SAEM(niter,nburnin,M,perdrix,paraminit,paramfix,delta)
    
    mul <- estim$mul[niter+1]
    muH <- estim$muH[niter+1]
    sigma2H <- estim$sigma2H[niter+1]
    sigma2l <- estim$sigma2l[niter+1]
    sigmaRes <- estim$sigmaRes[niter+1]
    
    biais[rep,] <- c(mul-mu_lambda,muH-mu_ch,sigma2l-sd2_lambda,sigma2H-sd2_ch,sigmaRes^2-sigma^2)
    std[rep,] <- sqrt(diag(solve(estim$FIM)))
  }
  
  res <- list(true.param=true.param,biais=biais,stdFIM=std)
  
  return(res)
}
