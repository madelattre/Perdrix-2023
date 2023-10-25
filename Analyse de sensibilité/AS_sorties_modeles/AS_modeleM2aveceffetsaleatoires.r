require(grDevices)
library(multisensi)

################################################################
## Modèle en 2D avec vigilance (cf articles bBaker et al, Ecological Modelling 2010; Billiard et al, JRSI 2018)

CE <- (sqrt(2)+log(1+sqrt(2)))/3
CV <- 2/3

esp2DM2<-function(D,ps,pv,cv,H,lambda){
    return((ps*(1-pv)*(sqrt(D)-1))/(pv*cv*(sqrt(D)-1)+CE*lambda*(1-pv)+H*ps*(1-pv)*(sqrt(D)-1)))
}

varmod2DM2<-function(D,ps,pv,cv,H,lambda,delta){
    ((1/delta)*ps*(1-pv)*(sqrt(D)-1)*(pv*ps*cv^2*(sqrt(D)-1)^2+(CV-CE^2)*lambda^2*ps*(1-pv)^2+(1-ps)*(pv*cv*(sqrt(D)-1)+CE*lambda*(1-pv))^2)/(pv*cv*(sqrt(D)-1)+CE*lambda*(1-pv)+H*ps*(1-pv)*(sqrt(D)-1))^3)
}


################################################################
################################################################
#Construction du plan d'expériences pour la moyenne, l'écart-type et le ratio espérance / écart-type de la réponse fonctionnelle

ASdesign=expand.grid(ps=seq(0.1,0.9,by=0.1), pv=seq(0.1,0.9,by=0.1), cv=c(0.3,1,5,10), H=c(0.5,2,5,10), lambda=c(0.1,1,5,10,100))
nbrun=dim(ASdesign)[1]

D <- c(1.1,2,3,5,10,15,25,50,100,200,400)
N <- 100 # nombre d'individus (perdrix)

ps_p=log(ASdesign$ps/(1-ASdesign$ps))

pv_p=log(ASdesign$pv/(1-ASdesign$pv))

mu_cv <- ASdesign$cv
sd_cv <- 0.1*mu_cv
logmu_cv=log(mu_cv)-log(1+sd_cv^2/((mu_cv)^2))/2
logsd_cv=sqrt(log(1+sd_cv^2/((mu_cv)^2)))

mu_H <- ASdesign$H
sd_H <- 0.1*mu_H
logmu_H=log(mu_H)-log(1+sd_H^2/((mu_H)^2))/2
logsd_H=sqrt(log(1+sd_H^2/((mu_H)^2)))

mu_lambda <- ASdesign$lambda
sd_lambda <- 0.1*mu_lambda
logmu_lambda=log(mu_lambda)-log(1+sd_lambda^2/((mu_lambda)^2))/2
logsd_lambda=sqrt(log(1+sd_lambda^2/((mu_lambda)^2)))

delta=10

y2DM2_m<-matrix(0, nrow = nbrun, ncol = length(D))
se2DM2_m<-matrix(0, nrow = nbrun, ncol = length(D))
noise2DM2_m<-matrix(0, nrow = nbrun, ncol = length(D))

#Simulations pour la matrice de plan d'expériences 
for (i in 1:nbrun)
{
    ps_AS <- rnorm(N,mean=ps_p[i],sd=0.1*abs(ps_p[i]))
    ps_AS <- exp(ps_AS)/(1+exp(ps_AS))
    
    pv_AS <- rnorm(N,mean=pv_p[i],sd=0.1*abs(pv_p[i]))
    pv_AS <- exp(pv_AS)/(1+exp(pv_AS))

    cv_AS <- exp(rnorm(N,mean=logmu_cv[i],sd=logsd_cv[i]))

    H_AS <- exp(rnorm(N,mean=logmu_H[i],sd=logsd_H[i]))
    
    lambda_AS <- exp(rnorm(N,mean=logmu_lambda[i],sd=logsd_lambda[i]))
    
    y2DM2=matrix(0, nrow = N, ncol = length(D))
    se2DM2=matrix(0, nrow = N, ncol = length(D))
    
    for (j in 1:N)
    {
    	y2DM2[j,]=esp2DM2(D,ps_AS[j],pv_AS[j],cv_AS[j],H_AS[j],lambda_AS[j])
        se2DM2[j,]=sqrt(varmod2DM2(D,ps_AS[j],pv_AS[j],cv_AS[j],H_AS[j],lambda_AS[j],delta))
    }
    y2DM2_m[i,] <- colMeans(y2DM2)
    se2DM2_m[i,] <- colMeans(se2DM2)
    noise2DM2_m[i,] <- y2DM2_m[i,]/se2DM2_m[i,]
}
################################################################

#Analyse de sensibilité des résultats - moyenne réponse fonctionnelle
repfonc_moy=data.frame(y2DM2_m)

ASrepfonc_moy.seq <- multisensi(model=repfonc_moy, design=ASdesign, reduction=NULL, center=FALSE)

plot(ASrepfonc_moy.seq, normalized = TRUE, color = terrain.colors, gsi.plot = FALSE)
title(main="Moyenne reponse fonctionnelle",xlab="D")

ASrepfonc_moy.pca <- multisensi(model=repfonc_moy, design=ASdesign, reduction=basis.ACP, scale=FALSE)
plot(ASrepfonc_moy.pca, graph = 1)
summary(ASrepfonc_moy.pca, digits = 2)
plot(ASrepfonc_moy.pca, graph = 2)
plot(ASrepfonc_moy.pca, graph = 3)

#---------------------------------------------------------------
#Analyse de sensibilité des résultats - écart-type réponse fonctionnelle
repfonc_std=data.frame(se2DM2_m)

ASrepfonc_std.seq <- multisensi(model=repfonc_std, design=ASdesign, reduction=NULL, center=FALSE)

plot(ASrepfonc_std.seq, normalized = TRUE, color = terrain.colors, gsi.plot = FALSE)
title(main="Ecart-type réponse fonctionnelle",xlab="D")

ASrepfonc_std.pca <- multisensi(model=repfonc_std, design=ASdesign, reduction=basis.ACP, scale=FALSE)
plot(ASrepfonc_std.pca, graph = 1)
summary(ASrepfonc_std.pca, digits = 2)
plot(ASrepfonc_std.pca, graph = 2)
plot(ASrepfonc_std.pca, graph = 3)

#---------------------------------------------------------------
#Analyse de sensibilité des résultats - ratio espérence / écart-type réponse fonctionnelle
repfonc_noise=data.frame(noise2DM2_m)

ASrepfonc_noise.seq <- multisensi(model=repfonc_noise, design=ASdesign, reduction=NULL, center=FALSE)

plot(ASrepfonc_noise.seq, normalized = TRUE, color = terrain.colors, gsi.plot = FALSE)
title(main="Ratio espérance/écart-type réponse fonctionnelle",xlab="D")

ASrepfonc_noise.pca <- multisensi(model=repfonc_noise, design=ASdesign, reduction=basis.ACP, scale=FALSE)
plot(ASrepfonc_noise.pca, graph = 1)
summary(ASrepfonc_noise.pca, digits = 2)
plot(ASrepfonc_noise.pca, graph = 2)
plot(ASrepfonc_noise.pca, graph = 3)


