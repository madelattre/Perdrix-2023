
require(grDevices)
library(multisensi)
library(fBasics)
library(matrixStats)
################################################################
#DEFINITION DU MODELE M1
################################################################

## Modèle M1 (sans vigilance) en 2D

CE <- (sqrt(2)+log(1+sqrt(2)))/3
CV <- 2/3

esp2DM1<-function(D,H,lambda){
  return(1/(CE*lambda*1/(sqrt(D)-1)+H))
}

varmod2DM1<-function(D,H,lambda,delta){
  return(1/delta*(CV-CE^2)*(lambda*1/(sqrt(D)-1))^2/(CE*lambda*1/(sqrt(D)-1)+H)^3)
}


################################################################
#ANALYSE DE SENSIBILITE MULTIVARIEE
################################################################

################################################################
#Construction du plan d'expériences pour les moyennes de lambda et H, ainsi que Delta. (Les écart-types de lambda et H sont déduits des moyennes.)
#Les outputs sont le moyenne, l'écart-type et le ratio espérance / écart-type de la réponse fonctionnelle.

ASdesign=expand.grid(H=c(0.5,1,2,5,10,50), lambda=c(0.1,1,5,10,100), delta=c(0.1,1,10,100))
nbrun=dim(ASdesign)[1]

D <- c(1.1,2,3,5,10,15,25,50,100,200,400)
N <- 100 # nombre d'individus (perdrix)

mu_H <- ASdesign$H
sd_H <- 0.1*mu_H
logmu_H=log(mu_H)-log(1+sd_H^2/((mu_H)^2))/2
logsd_H=sqrt(log(1+sd_H^2/((mu_H)^2)))

mu_lambda <- ASdesign$lambda
sd_lambda <- 0.1*mu_lambda
logmu_lambda=log(mu_lambda)-log(1+sd_lambda^2/((mu_lambda)^2))/2
logsd_lambda=sqrt(log(1+sd_lambda^2/((mu_lambda)^2)))

delta_AS=ASdesign$delta

y2DM1_m<-matrix(0, nrow = nbrun, ncol = length(D))
se2DM1_m<-matrix(0, nrow = nbrun, ncol = length(D))
noise2DM1_m<-matrix(0, nrow = nbrun, ncol = length(D))

#Simulations pour la matrice de plan d'expériences 
for (i in 1:nbrun)
{
    H_AS <- exp(rnorm(N,mean=logmu_H[i],sd=logsd_H[i]))
    
    lambda_AS <- exp(rnorm(N,mean=logmu_lambda[i],sd=logsd_lambda[i]))
    
    y2DM1=matrix(0, nrow = N, ncol = length(D))
    se2DM1=matrix(0, nrow = N, ncol = length(D))
    
    for (j in 1:N)
    {
    	y2DM1[j,]=esp2DM1(D,H_AS[j],lambda_AS[j])
        se2DM1[j,]=sqrt(varmod2DM1(D,H_AS[j],lambda_AS[j],delta_AS[i]))
    }
    y2DM1_m[i,] <- colMeans(y2DM1)
    se2DM1_m[i,] <- colMeans(se2DM1)
    noise2DM1_m[i,] <- y2DM1_m[i,]/se2DM1_m[i,]
}

################################################################
#Construction du plan d'expériences pour les moyennes et les écart-types de lambda et H, ainsi que Delta.
#Les outputs sont la moyenne, l'écart-type (et le ratio espérance / écart-type de la réponse fonctionnelle.)

#ASdesign=expand.grid(mu_H=c(0.5,1,2,5,10,50), sd_H=c(0.05,0.1,0.25,0.5), mu_lambda=c(0.1,1,5,10,100), sd_lambda=c(0.05,0.1,0.25,0.5), delta=c(0.1,1,10,100))

ASdesign=expand.grid(mu_H=c(0.5,1,2,5,10,50), sd_H=c(0.1,0.25,0.5,1,2,5), mu_lambda=c(0.1,1,5,10,100), sd_lambda=c(0.1,0.25,0.5,1,2,5), delta=c(0.1,1,10,100))

nbrun=dim(ASdesign)[1]

D <- c(1.1,1.5,2,2.5,3,4,5,6,7,8,9,10,15,25,50,100)
#D <- seq(1.1,3,0.1)

N <- 100 # nombre d'individus (perdrix)

mu_H <- ASdesign$mu_H
sd_H <- ASdesign$sd_H*mu_H
logmu_H=log(mu_H)-log(1+sd_H^2/((mu_H)^2))/2
logsd_H=sqrt(log(1+sd_H^2/((mu_H)^2)))

mu_lambda <- ASdesign$mu_lambda
sd_lambda <- ASdesign$sd_lambda*mu_lambda
logmu_lambda=log(mu_lambda)-log(1+sd_lambda^2/((mu_lambda)^2))/2
logsd_lambda=sqrt(log(1+sd_lambda^2/((mu_lambda)^2)))

delta_AS=ASdesign$delta

y2DM1_m<-matrix(0, nrow = nbrun, ncol = length(D))
se2DM1_m<-matrix(0, nrow = nbrun, ncol = length(D))
#noise2DM1_m<-matrix(0, nrow = nbrun, ncol = length(D))
y2DM1_v<-matrix(0, nrow = nbrun, ncol = length(D))
se2DM1_v<-matrix(0, nrow = nbrun, ncol = length(D))
#noise2DM1_v<-matrix(0, nrow = nbrun, ncol = length(D))


#Simulations pour la matrice de plan d'expériences
for (i in 1:nbrun) #boucle sur le nombre de répétitions
{
    H_AS <- exp(rnorm(N,mean=logmu_H[i],sd=logsd_H[i]))
    
    lambda_AS <- exp(rnorm(N,mean=logmu_lambda[i],sd=logsd_lambda[i]))
    
    y2DM1=matrix(0, nrow = N, ncol = length(D))
    se2DM1=matrix(0, nrow = N, ncol = length(D))
    
    for (j in 1:N) #boucle sur les individus
    {
        y2DM1[j,]=esp2DM1(D,H_AS[j],lambda_AS[j])
        se2DM1[j,]=sqrt(varmod2DM1(D,H_AS[j],lambda_AS[j],delta_AS[i]))
    }
    y2DM1_m[i,] <- colMeans(y2DM1)
    y2DM1_v[i,] <- colSds(y2DM1) #colStdevs(y2DM1)
    se2DM1_m[i,] <- colMeans(se2DM1)
    se2DM1_v[i,] <- colSds(se2DM1) #colStdevs(se2DM1)
    #    noise2DM1_m[i,] <- y2DM1_m[i,]/se2DM1_m[i,]
}

################################################################


#Analyse de sensibilité des résultats - moyenne réponse fonctionnelle - moyenne sur les effets aléatoires
repfonc_moy=data.frame(y2DM1_m)

ASrepfonc_moy.seq <- multisensi(model=repfonc_moy, design=ASdesign, reduction=NULL, center=FALSE)

plot(ASrepfonc_moy.seq, normalized = TRUE, color = terrain.colors, gsi.plot = FALSE)
title(main="M1 : Moyenne reponse fonctionnelle - moyenne",xlab="D")

ASrepfonc_moy.pca <- multisensi(model=repfonc_moy, design=ASdesign, reduction=basis.ACP, scale=FALSE)
plot(ASrepfonc_moy.pca, graph = 1)
summary(ASrepfonc_moy.pca, digits = 2)
plot(ASrepfonc_moy.pca, graph = 2)
plot(ASrepfonc_moy.pca, graph = 3)

#---------------------------------------------------------------

#Analyse de sensibilité des résultats - moyenne réponse fonctionnelle - écart-type sur les effets aléatoires

repfonc_sd=data.frame(y2DM1_v)

ASrepfonc_sd.seq <- multisensi(model=repfonc_sd, design=ASdesign, reduction=NULL, center=FALSE)

plot(ASrepfonc_sd.seq, normalized = TRUE, color = terrain.colors, gsi.plot = FALSE)
title(main="M1 : Moyenne reponse fonctionnelle - std",xlab="D")

ASrepfonc_sd.pca <- multisensi(model=repfonc_sd, design=ASdesign, reduction=basis.ACP, scale=FALSE)
plot(ASrepfonc_sd.pca, graph = 1)
summary(ASrepfonc_sd.pca, digits = 2)
plot(ASrepfonc_sd.pca, graph = 2)
plot(ASrepfonc_sd.pca, graph = 3)

#---------------------------------------------------------------

#Analyse de sensibilité des résultats - écart-type réponse fonctionnelle - moyenne sur les effets aléatoires
repfonc_std_moy=data.frame(se2DM1_m)

ASrepfonc_std_moy.seq <- multisensi(model=repfonc_std_moy, design=ASdesign, reduction=NULL, center=FALSE)

plot(ASrepfonc_std_moy.seq, normalized = TRUE, color = terrain.colors, gsi.plot = FALSE)
title(main="M1 : Ecart-type réponse fonctionnelle - moyenne",xlab="D")

ASrepfonc_std_moy.pca <- multisensi(model=repfonc_std_moy, design=ASdesign, reduction=basis.ACP, scale=FALSE)
plot(ASrepfonc_std_moy.pca, graph = 1)
summary(ASrepfonc_std_moy.pca, digits = 2)
plot(ASrepfonc_std_moy.pca, graph = 2)
plot(ASrepfonc_std_moy.pca, graph = 3)

#---------------------------------------------------------------

#Analyse de sensibilité des résultats - écart-type réponse fonctionnelle - écart-type sur les effets aléatoires
repfonc_std_std=data.frame(se2DM1_v)

ASrepfonc_std_std.seq <- multisensi(model=repfonc_std_std, design=ASdesign, reduction=NULL, center=FALSE)

plot(ASrepfonc_std_std.seq, normalized = TRUE, color = terrain.colors, gsi.plot = FALSE)
title(main="M1 : Ecart-type réponse fonctionnelle - std",xlab="D")

ASrepfonc_std_std.pca <- multisensi(model=repfonc_std_std, design=ASdesign, reduction=basis.ACP, scale=FALSE)
plot(ASrepfonc_std_std.pca, graph = 1)
summary(ASrepfonc_std_std.pca, digits = 2)
plot(ASrepfonc_std_std.pca, graph = 2)
plot(ASrepfonc_std_std.pca, graph = 3)


################################################################
#Construction du plan d'expériences pour les moyennes et les écart-types de lambda et H, ainsi que Delta.
#Les outputs sont le ratio espérance / écart-type de la réponse fonctionnelle (moyenne des ratios ou ratio des moyennes sur tous les individus).

ASdesign=expand.grid(mu_H=c(0.5,1,2,5,10,50), sd_H=c(0.1,0.25,0.5,1,2,5), mu_lambda=c(0.1,1,5,10,100), sd_lambda=c(0.1,0.25,0.5,1,2,5), delta=c(0.1,1,10,100))

nbrun=dim(ASdesign)[1]

D <- c(1.1,1.5,2,2.5,3,4,5,6,7,8,9,10,15,25,50,100)
#D <- seq(1.1,3,0.1)

N <- 100 # nombre d'individus (perdrix)

mu_H <- ASdesign$mu_H
sd_H <- ASdesign$sd_H*mu_H
logmu_H=log(mu_H)-log(1+sd_H^2/((mu_H)^2))/2
logsd_H=sqrt(log(1+sd_H^2/((mu_H)^2)))

mu_lambda <- ASdesign$mu_lambda
sd_lambda <- ASdesign$sd_lambda*mu_lambda
logmu_lambda=log(mu_lambda)-log(1+sd_lambda^2/((mu_lambda)^2))/2
logsd_lambda=sqrt(log(1+sd_lambda^2/((mu_lambda)^2)))

delta_AS=ASdesign$delta

y2DM1_m<-matrix(0, nrow = nbrun, ncol = length(D))
se2DM1_m<-matrix(0, nrow = nbrun, ncol = length(D))
noise2DM1_rm<-matrix(0, nrow = nbrun, ncol = length(D))
noise2DM1_mr<-matrix(0, nrow = nbrun, ncol = length(D))

#y2DM1_v<-matrix(0, nrow = nbrun, ncol = length(D))
#se2DM1_v<-matrix(0, nrow = nbrun, ncol = length(D))
#noise2DM1_v<-matrix(0, nrow = nbrun, ncol = length(D))


#Simulations pour la matrice de plan d'expériences
for (i in 1:nbrun) #boucle sur le nombre d'éléments dans le plan d'expérience
{
    H_AS <- exp(rnorm(N,mean=logmu_H[i],sd=logsd_H[i]))
    
    lambda_AS <- exp(rnorm(N,mean=logmu_lambda[i],sd=logsd_lambda[i]))
    
    y2DM1=matrix(0, nrow = N, ncol = length(D))
    se2DM1=matrix(0, nrow = N, ncol = length(D))
    
    for (j in 1:N) #boucle sur les individus
    {
        y2DM1[j,]=esp2DM1(D,H_AS[j],lambda_AS[j])
        se2DM1[j,]=sqrt(varmod2DM1(D,H_AS[j],lambda_AS[j],delta_AS[i]))
    }
    y2DM1_m[i,] <- colMeans(y2DM1)
    #y2DM1_v[i,] <- colStdevs(y2DM1)
    se2DM1_m[i,] <- colMeans(se2DM1)
    #se2DM1_v[i,] <- colStdevs(se2DM1)
    noise2DM1_rm[i,] <- y2DM1_m[i,]/se2DM1_m[i,]
    noise2DM1_mr[i,] <- colMeans(y2DM1/se2DM1)
}

################################################################

#Analyse de sensibilité des résultats - ratio signal sur bruit comme ratio entre la moyenne de la réponse fonctionnele moyennée sur tous les individus et l'écart-type de la réponse fonctionnelle, moyenné sur tous les individus
repfonc_moy=data.frame(noise2DM1_rm)

ASrepfonc_moy.seq <- multisensi(model=repfonc_moy, design=ASdesign, reduction=NULL, center=FALSE)

plot(ASrepfonc_moy.seq, normalized = TRUE, color = terrain.colors, gsi.plot = FALSE)
title(main="M1 : E[R]/std[R], num. et denomin. moyennés sur tous les individus",xlab="D")

ASrepfonc_moy.pca <- multisensi(model=repfonc_moy, design=ASdesign, reduction=basis.ACP, scale=FALSE)
plot(ASrepfonc_moy.pca, graph = 1)
summary(ASrepfonc_moy.pca, digits = 2)
plot(ASrepfonc_moy.pca, graph = 2)
plot(ASrepfonc_moy.pca, graph = 3)


#---------------------------------------------------------------
#Analyse de sensibilité des résultats - ratio signal sur bruit comme ratio entre la moyenne de la réponse fonctionnele et l'écart-type de la réponse fonctionnelle, ratio moyenné sur tous les individus
repfonc_moy=data.frame(noise2DM1_mr)

ASrepfonc_moy.seq <- multisensi(model=repfonc_moy, design=ASdesign, reduction=NULL, center=FALSE)

plot(ASrepfonc_moy.seq, normalized = TRUE, color = terrain.colors, gsi.plot = FALSE)
title(main="M1 : moyenne de E[R]_i/std[R]_i, i est l'individu",xlab="D")

ASrepfonc_moy.pca <- multisensi(model=repfonc_moy, design=ASdesign, reduction=basis.ACP, scale=FALSE)
plot(ASrepfonc_moy.pca, graph = 1)
summary(ASrepfonc_moy.pca, digits = 2)
plot(ASrepfonc_moy.pca, graph = 2)
plot(ASrepfonc_moy.pca, graph = 3)

################################################################


#---------------------------------------------------------------
#D <- c(1.1,2,3,5,10,15,25,50,100,200,400)

#lambda <- 8
#H <- 0.5
#delta=10

D <- c(1.1,1.2,1.3,1.5,2,3,5,10,15,25)

lambda <- 5
H <- 5
delta=1

y2DM1=esp2DM1(D,H,lambda)
se2DM1=sqrt(varmod2DM1(D,H,lambda,delta))
plot(y2DM1~D,type="o",ylim=c(0,1.5))
arrows(D,y2DM1,D,y2DM1+1.96*se2DM1,length=0.05,angle=90)
arrows(D,y2DM1,D,y2DM1-1.96*se2DM1,length=0.05,angle=90)

esp_lambda=lambda; var_lambda=(0.1*esp_lambda)^2
#var_lambda=(esp_lambda)^2
esp_logl=log(esp_lambda)-0.5*log(1+var_lambda/((esp_lambda)^2))
sd_logl=sqrt(log(1+var_lambda/((esp_lambda)^2)))

esp_H=H; var_H=(0.1*H)^2
#var_H=(esp_H)^2
esp_logH=log(esp_H)-0.5*log(1+var_H/((esp_H)^2))
sd_logH=sqrt(log(1+var_H/((esp_H)^2)))


for (i in 2:5)
{
lambda=exp(rnorm(1,esp_logl,sd_logl))
H=exp(rnorm(1,esp_logH,sd_logH))
par(new=T)
y2DM1=esp2DM1(D,H,lambda)
se2DM1=sqrt(varmod2DM1(D,H,lambda,delta))
lines(y2DM1~D,type="o",pch=i)
arrows(D,y2DM1,D,y2DM1+1.96*se2DM1,length=0.05,angle=90)
arrows(D,y2DM1,D,y2DM1-1.96*se2DM1,length=0.05,angle=90)
}

################################################################
#ANALYSE DE SENSIBILITE UNIVARIEE
################################################################


#AS univariée mu_H
#Variante avec les écart-types de lambda et H
#ASdesign=expand.grid(mu_H=c(0.5,1,2,5,10,50), sd_H=c(0.05,0.1,0.25,0.5), mu_lambda=c(0.1,1,5,10,100), sd_lambda=c(0.05,0.1,0.25,0.5), delta=c(0.1,1,10,100))

mu_H=c(0.5,1,2,5,10,50)
#sd_H=c(0.1,0.25,0.5,1,2,5)
#mu_lambda=c(0.1,1,5,10,100)
#sd_lambda=c(0.1,0.25,0.5,1,2,5)
#delta=c(0.1,1,10,100))

nbrun=length(mu_H)

D <- c(1.1,1.5,2,2.5,3,4,5,6,7,8,9,10,15,25,50,100)
#D <- seq(1.1,3,0.1)

N <- 100 # nombre d'individus (perdrix)

sd_H <- 0.25
logmu_H=log(mu_H)-log(1+sd_H^2/((mu_H)^2))/2
logsd_H=sqrt(log(1+sd_H^2/((mu_H)^2)))

mu_lambda <- 8
sd_lambda <- 0.25
logmu_lambda=log(mu_lambda)-log(1+sd_lambda^2/((mu_lambda)^2))/2
logsd_lambda=sqrt(log(1+sd_lambda^2/((mu_lambda)^2)))

delta_AS=10

y2DM1_m<-matrix(0, nrow = nbrun, ncol = length(D))
se2DM1_m<-matrix(0, nrow = nbrun, ncol = length(D))
#noise2DM1_m<-matrix(0, nrow = nbrun, ncol = length(D))
y2DM1_v<-matrix(0, nrow = nbrun, ncol = length(D))
se2DM1_v<-matrix(0, nrow = nbrun, ncol = length(D))
#noise2DM1_v<-matrix(0, nrow = nbrun, ncol = length(D))


#Simulations pour la matrice de plan d'expériences
for (i in 1:nbrun)
{
    H_AS <- exp(rnorm(N,mean=logmu_H[i],sd=logsd_H))
    
    lambda_AS <- exp(rnorm(N,mean=logmu_lambda,sd=logsd_lambda))
    
    y2DM1=matrix(0, nrow = N, ncol = length(D))
    se2DM1=matrix(0, nrow = N, ncol = length(D))
    
    for (j in 1:N)
    {
        y2DM1[j,]=esp2DM1(D,H_AS[j],lambda_AS[j])
        se2DM1[j,]=sqrt(varmod2DM1(D,H_AS[j],lambda_AS[j],delta_AS))
    }
    y2DM1_m[i,] <- colMeans(y2DM1)
    y2DM1_v[i,] <- colStdevs(y2DM1)
    se2DM1_m[i,] <- colMeans(se2DM1)
    se2DM1_v[i,] <- colStdevs(se2DM1)
    #    noise2DM1_m[i,] <- y2DM1_m[i,]/se2DM1_m[i,]
}

#Plot Reponse fonctionnelle moy / écart-type du processus
plot(y2DM1_m[1,]~D,type="o",pch=1,ylim=c(0,1),col=colors()[1*10],main="H")
arrows(D,y2DM1_m[1,],D,y2DM1_m[1,]+1.96*se2DM1_m[1,],length=0.05,angle=90,col=colors()[1*10])
arrows(D,y2DM1_m[1,],D,y2DM1_m[1,]-1.96*se2DM1_m[1,],length=0.05,angle=90,col=colors()[1*10])

for (i in 2:nbrun)
{
    lines(y2DM1_m[i,]~D,type="o",col=colors()[i*10],pch=i)
    arrows(D,y2DM1_m[i,],D,y2DM1_m[i,]+1.96*se2DM1_m[i,],length=0.05,angle=90,col=colors()[i*10])
    arrows(D,y2DM1_m[i,],D,y2DM1_m[i,]-1.96*se2DM1_m[i,],length=0.05,angle=90,col=colors()[i*10])
}

#Plot Reponse fonctionnelle moy / écart-type sur les individus
plot(y2DM1_m[1,]~D,type="o",pch=1,ylim=c(0,1.5),col=colors()[1*10],main="H")
arrows(D,y2DM1_m[1,],D,y2DM1_m[1,]+1.96*sqrt(y2DM1_v[1,]),length=0.05,angle=90,col=colors()[1*10])
arrows(D,y2DM1_m[1,],D,y2DM1_m[1,]-1.96*sqrt(y2DM1_v[1,]),length=0.05,angle=90,col=colors()[1*10])

for (i in 2:nbrun)
{
    lines(y2DM1_m[i,]~D,type="o",col=colors()[i*10],pch=i)
    arrows(D,y2DM1_m[i,],D,y2DM1_m[i,]+1.96*sqrt(y2DM1_v[i,]),length=0.05,angle=90,col=colors()[i*10])
    arrows(D,y2DM1_m[i,],D,y2DM1_m[i,]-1.96*sqrt(y2DM1_v[i,]),length=0.05,angle=90,col=colors()[i*10])
}

#Plot Reponse fonctionnelle écart-type / écart-type sur les individus
plot(se2DM1_m[1,]~D,type="o",pch=1,ylim=c(0,0.3),col=colors()[1*10],main="H")
arrows(D,se2DM1_m[1,],D,se2DM1_m[1,]+1.96*sqrt(se2DM1_v[1,]),length=0.05,angle=90,col=colors()[1*10])
arrows(D,se2DM1_m[1,],D,se2DM1_m[1,]-1.96*sqrt(se2DM1_v[1,]),length=0.05,angle=90,col=colors()[1*10])

for (i in 2:nbrun)
{
    lines(se2DM1_m[i,]~D,type="o",col=colors()[i*10],pch=i)
    arrows(D,se2DM1_m[i,],D,se2DM1_m[i,]+1.96*sqrt(se2DM1_v[i,]),length=0.05,angle=90,col=colors()[i*10])
    arrows(D,se2DM1_m[i,],D,se2DM1_m[i,]-1.96*sqrt(se2DM1_v[i,]),length=0.05,angle=90,col=colors()[i*10])
}


################################################################
#AS univariée mu_lambda
#Variante avec les écart-types de lambda et H
#ASdesign=expand.grid(mu_H=c(0.5,1,2,5,10,50), sd_H=c(0.05,0.1,0.25,0.5), mu_lambda=c(0.1,1,5,10,100), sd_lambda=c(0.05,0.1,0.25,0.5), delta=c(0.1,1,10,100))

#mu_H=c(0.5,1,2,5,10,50)
#sd_H=c(0.1,0.25,0.5,1,2,5)
mu_lambda=c(0.1,1,5,10,100)
#sd_lambda=c(0.1,0.25,0.5,1,2,5)
#delta=c(0.1,1,10,100))

nbrun=length(mu_lambda)

D <- c(1.1,1.5,2,2.5,3,4,5,6,7,8,9,10,15,25,50,100)
#D <- seq(1.1,3,0.1)

N <- 100 # nombre d'individus (perdrix)

mu_H=0.5
sd_H <- 0.25
logmu_H=log(mu_H)-log(1+sd_H^2/((mu_H)^2))/2
logsd_H=sqrt(log(1+sd_H^2/((mu_H)^2)))

sd_lambda <- 0.25
logmu_lambda=log(mu_lambda)-log(1+sd_lambda^2/((mu_lambda)^2))/2
logsd_lambda=sqrt(log(1+sd_lambda^2/((mu_lambda)^2)))

delta_AS=10

y2DM1_m<-matrix(0, nrow = nbrun, ncol = length(D))
se2DM1_m<-matrix(0, nrow = nbrun, ncol = length(D))
#noise2DM1_m<-matrix(0, nrow = nbrun, ncol = length(D))
y2DM1_v<-matrix(0, nrow = nbrun, ncol = length(D))
se2DM1_v<-matrix(0, nrow = nbrun, ncol = length(D))
#noise2DM1_v<-matrix(0, nrow = nbrun, ncol = length(D))


#Simulations pour la matrice de plan d'expériences
for (i in 1:nbrun)
{
    H_AS <- exp(rnorm(N,mean=logmu_H,sd=logsd_H))
    
    lambda_AS <- exp(rnorm(N,mean=logmu_lambda[i],sd=logsd_lambda))
    
    y2DM1=matrix(0, nrow = N, ncol = length(D))
    se2DM1=matrix(0, nrow = N, ncol = length(D))
    
    for (j in 1:N)
    {
        y2DM1[j,]=esp2DM1(D,H_AS[j],lambda_AS[j])
        se2DM1[j,]=sqrt(varmod2DM1(D,H_AS[j],lambda_AS[j],delta_AS))
    }
    y2DM1_m[i,] <- colMeans(y2DM1)
    y2DM1_v[i,] <- colStdevs(y2DM1)
    se2DM1_m[i,] <- colMeans(se2DM1)
    se2DM1_v[i,] <- colStdevs(se2DM1)
    #    noise2DM1_m[i,] <- y2DM1_m[i,]/se2DM1_m[i,]
}

#Plot Reponse fonctionnelle moy / écart-type du processus
plot(y2DM1_m[1,]~D,type="o",pch=1,ylim=c(0,3),col=colors()[1*10],main="lambda")
arrows(D,y2DM1_m[1,],D,y2DM1_m[1,]+1.96*se2DM1_m[1,],length=0.05,angle=90,col=colors()[1*10])
arrows(D,y2DM1_m[1,],D,y2DM1_m[1,]-1.96*se2DM1_m[1,],length=0.05,angle=90,col=colors()[1*10])

for (i in 2:nbrun)
{
    lines(y2DM1_m[i,]~D,type="o",col=colors()[i*10],pch=i)
    arrows(D,y2DM1_m[i,],D,y2DM1_m[i,]+1.96*se2DM1_m[i,],length=0.05,angle=90,col=colors()[i*10])
    arrows(D,y2DM1_m[i,],D,y2DM1_m[i,]-1.96*se2DM1_m[i,],length=0.05,angle=90,col=colors()[i*10])
}

#Plot Reponse fonctionnelle moy / écart-type sur les individus
plot(y2DM1_m[1,]~D,type="o",pch=1,ylim=c(0,5),col=colors()[1*10],main="lambda")
arrows(D,y2DM1_m[1,],D,y2DM1_m[1,]+1.96*sqrt(y2DM1_v[1,]),length=0.05,angle=90,col=colors()[1*10])
arrows(D,y2DM1_m[1,],D,y2DM1_m[1,]-1.96*sqrt(y2DM1_v[1,]),length=0.05,angle=90,col=colors()[1*10])

for (i in 2:nbrun)
{
    lines(y2DM1_m[i,]~D,type="o",col=colors()[i*10],pch=i)
    arrows(D,y2DM1_m[i,],D,y2DM1_m[i,]+1.96*sqrt(y2DM1_v[i,]),length=0.05,angle=90,col=colors()[i*10])
    arrows(D,y2DM1_m[i,],D,y2DM1_m[i,]-1.96*sqrt(y2DM1_v[i,]),length=0.05,angle=90,col=colors()[i*10])
}

#Plot Reponse fonctionnelle variance / écart-type sur les individus
plot(se2DM1_m[1,]~D,type="o",pch=1,ylim=c(0,0.5),col=colors()[1*10],main="lambda")
arrows(D,se2DM1_m[1,],D,se2DM1_m[1,]+1.96*sqrt(se2DM1_v[1,]),length=0.05,angle=90,col=colors()[1*10])
arrows(D,se2DM1_m[1,],D,se2DM1_m[1,]-1.96*sqrt(se2DM1_v[1,]),length=0.05,angle=90,col=colors()[1*10])

for (i in 2:nbrun)
{
    lines(se2DM1_m[i,]~D,type="o",col=colors()[i*10],pch=i)
    arrows(D,se2DM1_m[i,],D,se2DM1_m[i,]+1.96*sqrt(se2DM1_v[i,]),length=0.05,angle=90,col=colors()[i*10])
    arrows(D,se2DM1_m[i,],D,se2DM1_m[i,]-1.96*sqrt(se2DM1_v[i,]),length=0.05,angle=90,col=colors()[i*10])
}

#------------------------------------------
#Analyse univariée lambda modèle à effets fixes
#------------------------------------------

lambda=c(0.1,1,5,10,100)
H=0.5
delta=1

nbrun=length(lambda)

D <- c(1.1,1.5,2,2.5,3,4,5,6,7,8,9,10,15,25,50,100)

y2DM1_f<-matrix(0, nrow = nbrun, ncol = length(D))
se2DM1_f<-matrix(0, nrow = nbrun, ncol = length(D))

for (i in 1:nbrun)
{
    y2DM1_f[i,]=esp2DM1(D,H,lambda[i])
   se2DM1_f[i,]=sqrt(varmod2DM1(D,H,lambda[i],delta))
}

x11()
plot(y2DM1_f[1,]~D,type="o",pch=1,ylim=c(0,2),col=colors()[1*10],main="Analyse sensibilité lambda",ylab="Moyenne rép. fonc.")
arrows(D,y2DM1_f[1,],D,y2DM1_f[1,]+1.96*se2DM1_f[1,],length=0.05,angle=90,col=colors()[1*10])
arrows(D,y2DM1_f[1,],D,y2DM1_f[1,]-1.96*se2DM1_f[1,],length=0.05,angle=90,col=colors()[1*10])

for (i in 2:nbrun)
{
    lines(y2DM1_f[i,]~D,type="o",pch=1,ylim=c(0,2),col=colors()[i*10])
    arrows(D,y2DM1_f[i,],D,y2DM1_f[i,]+1.96*se2DM1_f[i,],length=0.05,angle=90,col=colors()[i*10])
    arrows(D,y2DM1_f[i,],D,y2DM1_f[i,]-1.96*se2DM1_f[i,],length=0.05,angle=90,col=colors()[i*10])
}

x11()
plot(se2DM1_f[1,]~D,type="o",pch=1,ylim=c(0,0.2),col=colors()[1*10],main="Analyse sensibilité lambda",ylab="Ecart-type rép. fonc.")

for (i in 2:nbrun)
{
    lines(se2DM1_f[i,]~D,type="o",pch=1,ylim=c(0,0.2),col=colors()[i*10])
  }

#------------------------------------------
#Analyse univariée H modèle à effets fixes
#------------------------------------------
H=c(0.5,1,2,5,10,50)
lambda=8
delta=1

nbrun=length(H)

D <- c(1.1,1.5,2,2.5,3,4,5,6,7,8,9,10,15,25,50,100)

y2DM1_f<-matrix(0, nrow = nbrun, ncol = length(D))
se2DM1_f<-matrix(0, nrow = nbrun, ncol = length(D))

for (i in 1:nbrun)
{
    y2DM1_f[i,]=esp2DM1(D,H[i],lambda)
    se2DM1_f[i,]=sqrt(varmod2DM1(D,H[i],lambda,delta))
}

x11()
plot(y2DM1_f[1,]~D,type="o",pch=1,ylim=c(0,1.5),col=colors()[1*10],main="Analyse sensibilité H",ylab="Moyenne rép. fonc.")
arrows(D,y2DM1_f[1,],D,y2DM1_f[1,]+1.96*se2DM1_f[1,],length=0.05,angle=90,col=colors()[1*10])
arrows(D,y2DM1_f[1,],D,y2DM1_f[1,]-1.96*se2DM1_f[1,],length=0.05,angle=90,col=colors()[1*10])

for (i in 2:nbrun)
{
    lines(y2DM1_f[i,]~D,type="o",pch=1,ylim=c(0,1.5),col=colors()[i*10])
    arrows(D,y2DM1_f[i,],D,y2DM1_f[i,]+1.96*se2DM1_f[i,],length=0.05,angle=90,col=colors()[i*10])
    arrows(D,y2DM1_f[i,],D,y2DM1_f[i,]-1.96*se2DM1_f[i,],length=0.05,angle=90,col=colors()[i*10])
}

x11()
plot(se2DM1_f[1,]~D,type="o",pch=1,ylim=c(0,0.2),col=colors()[1*10],main="Analyse sensibilité H",ylab="Ecart-type rép. fonc.")

for (i in 2:nbrun)
{
    lines(se2DM1_f[i,]~D,type="o",pch=1,ylim=c(0,0.2),col=colors()[i*10])
}


################################################################
#AS univariée Delta
#Variante avec les écart-types de lambda et H
#ASdesign=expand.grid(mu_H=c(0.5,1,2,5,10,50), sd_H=c(0.05,0.1,0.25,0.5), mu_lambda=c(0.1,1,5,10,100), sd_lambda=c(0.05,0.1,0.25,0.5), delta=c(0.1,1,10,100))

#mu_H=c(0.5,1,2,5,10,50)
#sd_H=c(0.1,0.25,0.5,1,2,5)
#mu_lambda=c(0.1,1,5,10,100)
#sd_lambda=c(0.1,0.25,0.5,1,2,5)
delta=c(0.1,1,10,100)

nbrun=length(delta)

D <- c(1.1,1.5,2,2.5,3,4,5,6,7,8,9,10,15,25,50,100)
#D <- seq(1.1,3,0.1)

N <- 100 # nombre d'individus (perdrix)


mu_H=0.5
sd_H <- 0.25
logmu_H=log(mu_H)-log(1+sd_H^2/((mu_H)^2))/2
logsd_H=sqrt(log(1+sd_H^2/((mu_H)^2)))

mu_lambda=8
sd_lambda <- 0.25
logmu_lambda=log(mu_lambda)-log(1+sd_lambda^2/((mu_lambda)^2))/2
logsd_lambda=sqrt(log(1+sd_lambda^2/((mu_lambda)^2)))

delta_AS=delta

y2DM1_m<-matrix(0, nrow = nbrun, ncol = length(D))
se2DM1_m<-matrix(0, nrow = nbrun, ncol = length(D))
#noise2DM1_m<-matrix(0, nrow = nbrun, ncol = length(D))
y2DM1_v<-matrix(0, nrow = nbrun, ncol = length(D))
se2DM1_v<-matrix(0, nrow = nbrun, ncol = length(D))
#noise2DM1_v<-matrix(0, nrow = nbrun, ncol = length(D))


#Simulations pour la matrice de plan d'expériences
for (i in 1:nbrun)
{
    H_AS <- exp(rnorm(N,mean=logmu_H,sd=logsd_H))
    
    lambda_AS <- exp(rnorm(N,mean=logmu_lambda,sd=logsd_lambda))
    
    y2DM1=matrix(0, nrow = N, ncol = length(D))
    se2DM1=matrix(0, nrow = N, ncol = length(D))
    
    for (j in 1:N)
    {
        y2DM1[j,]=esp2DM1(D,H_AS[j],lambda_AS[j])
        se2DM1[j,]=sqrt(varmod2DM1(D,H_AS[j],lambda_AS[j],delta_AS[i]))
    }
    y2DM1_m[i,] <- colMeans(y2DM1)
    y2DM1_v[i,] <- colStdevs(y2DM1)
    se2DM1_m[i,] <- colMeans(se2DM1)
    se2DM1_v[i,] <- colStdevs(se2DM1)
    #    noise2DM1_m[i,] <- y2DM1_m[i,]/se2DM1_m[i,]
}

#Plot Reponse fonctionnelle moy / écart-type du processus
plot(y2DM1_m[1,]~D,type="o",pch=1,ylim=c(0,2),col=colors()[1*10],main="Delta")
arrows(D,y2DM1_m[1,],D,y2DM1_m[1,]+1.96*se2DM1_m[1,],length=0.05,angle=90,col=colors()[1*10])
arrows(D,y2DM1_m[1,],D,y2DM1_m[1,]-1.96*se2DM1_m[1,],length=0.05,angle=90,col=colors()[1*10])

for (i in 2:nbrun)
{
    lines(y2DM1_m[i,]~D,type="o",col=colors()[i*10],pch=i)
    arrows(D,y2DM1_m[i,],D,y2DM1_m[i,]+1.96*se2DM1_m[i,],length=0.05,angle=90,col=colors()[i*10])
    arrows(D,y2DM1_m[i,],D,y2DM1_m[i,]-1.96*se2DM1_m[i,],length=0.05,angle=90,col=colors()[i*10])
}

#Plot Reponse fonctionnelle moy / écart-type sur les individus
plot(y2DM1_m[1,]~D,type="o",pch=1,ylim=c(0,2),col=colors()[1*10],main="Delta")
arrows(D,y2DM1_m[1,],D,y2DM1_m[1,]+1.96*sqrt(y2DM1_v[1,]),length=0.05,angle=90,col=colors()[1*10])
arrows(D,y2DM1_m[1,],D,y2DM1_m[1,]-1.96*sqrt(y2DM1_v[1,]),length=0.05,angle=90,col=colors()[1*10])

for (i in 2:nbrun)
{
    lines(y2DM1_m[i,]~D,type="o",col=colors()[i*10],pch=i)
    arrows(D,y2DM1_m[i,],D,y2DM1_m[i,]+1.96*sqrt(y2DM1_v[i,]),length=0.05,angle=90,col=colors()[i*10])
    arrows(D,y2DM1_m[i,],D,y2DM1_m[i,]-1.96*sqrt(y2DM1_v[i,]),length=0.05,angle=90,col=colors()[i*10])
}

#Plot Reponse fonctionnelle variance / écart-type sur les individus
plot(se2DM1_m[1,]~D,type="o",pch=1,ylim=c(0,1.5),col=colors()[1*10],main="Delta")
arrows(D,se2DM1_m[1,],D,se2DM1_m[1,]+1.96*sqrt(se2DM1_v[1,]),length=0.05,angle=90,col=colors()[1*10])
arrows(D,se2DM1_m[1,],D,se2DM1_m[1,]-1.96*sqrt(se2DM1_v[1,]),length=0.05,angle=90,col=colors()[1*10])

for (i in 2:nbrun)
{
    lines(se2DM1_m[i,]~D,type="o",col=colors()[i*10],pch=i)
    arrows(D,se2DM1_m[i,],D,se2DM1_m[i,]+1.96*sqrt(se2DM1_v[i,]),length=0.05,angle=90,col=colors()[i*10])
    arrows(D,se2DM1_m[i,],D,se2DM1_m[i,]-1.96*sqrt(se2DM1_v[i,]),length=0.05,angle=90,col=colors()[i*10])
}



###################
#Scenario 3
#mu_H = 5; std_H=0.25*mu_H; mu_lambda = 1; std_lambda = 0.25 * mu_lambda ; delta = 1 ;
#Scenario 4
#mu_H = 5; std_H=0.1*mu_H; mu_lambda = 5; std_lambda = 0.1 * mu_lambda ; Delta = 1;
#Scenario 1
mu_H = 0.5 ; std_H=0.25*mu_H ; mu_lambda = 8 ; std_lambda = 0.25 * mu_lambda ; delta = 10 ;
#Scenario X
mu_H = 0.5 ; std_H=0.75*mu_H ; mu_lambda = 1 ; std_lambda = 0.25 * mu_lambda ; delta = 1 ;
mu_H = 10 ; std_H=0.5*mu_H ; mu_lambda = 20 ; std_lambda = mu_lambda ; delta = 1 ;
mu_H = 0.5 ; sd_H=0.25 ; mu_lambda = 1 ; sd_lambda = 0.5 ; delta = 1 ;

#D <- c(1.1,1.5,2,2.5,3,4,5,6,7,8,9,10,15,25,50,100)
#Scenarios 3 & 4
#D <-c(1.1,1.2,1.3,1.5,2,3,5,10,15,25)
#Scenarios 1 & 2
D <- c(1.1,2,3,5,10,15,25,50,100,200,400)
D <- c(1.1,2,3,5,10,15,25,50,100,200,400)
D <- c(1.1,2,3,5,7,10,15,25,50,100,200)


N <- 100 # nombre d'individus (perdrix)

logmu_H=log(mu_H)-log(1+sd_H^2/((mu_H)^2))/2
logsd_H=sqrt(log(1+sd_H^2/((mu_H)^2)))

logmu_lambda=log(mu_lambda)-log(1+sd_lambda^2/((mu_lambda)^2))/2
logsd_lambda=sqrt(log(1+sd_lambda^2/((mu_lambda)^2)))

H_AS <- exp(rnorm(N,mean=logmu_H,sd=logsd_H))
lambda_AS <- exp(rnorm(N,mean=logmu_lambda,sd=logsd_lambda))

y2DM1<-matrix(0, nrow = N, ncol = length(D))
se2DM1<-matrix(0, nrow = N, ncol = length(D))

for (j in 1:N)
{
    y2DM1[j,]=esp2DM1(D,H_AS[j],lambda_AS[j])
    se2DM1[j,]=sqrt(varmod2DM1(D,H_AS[j],lambda_AS[j],delta))
}

y2DM1_m <- colMeans(y2DM1)
y2DM1_v <- colStdevs(y2DM1)
se2DM1_m <- colMeans(se2DM1)
se2DM1_v <- colStdevs(se2DM1)


plot(y2DM1_m~D,type="o",pch=1,ylim=c(0,2.5),col=colors()[1*10],main="Moyenne Réponse fonctionnelle")
arrows(D,y2DM1_m,D,y2DM1_m+1.96*se2DM1_m,length=0.05,angle=90,col=colors()[1*10])
arrows(D,y2DM1_m,D,y2DM1_m-1.96*se2DM1_m,length=0.05,angle=90,col=colors()[1*10])

