require(grDevices)

################################################################
## Modèle en 2D

CE <- (sqrt(2)+log(1+sqrt(2)))/3
CV <- 2/3

esp2DM2<-function(D,ps,pv,cv,H,lambda){
  return((ps*(1-pv)*(sqrt(D)-1))/(pv*CV*(sqrt(D)-1)+CE*lambda*(1-pv)+H*ps*(1-pv)*(sqrt(D)-1)))
}
#Le cv n'intervient pas, on peux l'enlever des arguments de la fonction esp2DM2

varmod2DM2<-function(D,ps,pv,cv,H,lambda,delta){
  ((1/delta)*ps*(1-pv)*(sqrt(D)-1)*(pv*ps*cv^2*(sqrt(D)-1)^2+(CV-CE^2)*lambda^2*ps*(1-pv)^2+(1-ps)*(pv*cv*(sqrt(D)-1)+CE*lambda*(1-pv))^2)/(pv*cv*(sqrt(D)-1)+CE*lambda*(1-pv)+H*ps*(1-pv)*(sqrt(D)-1))^3)
}

################################################################
################################################################
####Plot Modele 2DM2
#D <- c(1.5,2,3,5,10,25,50,75,100)
#D <- c(1.5,2,3,5,10,25,50)
#delta=60*4.5

#D = c(5,10,15,25,50,100,200,400) #Valeurs article Baker et al, 2010
D = c(1.5,2,3,5,10,100,200,400)
#D <- c(1.5,2,3,5,10)
ps=1
pv=0.5
cv=0.37
H=7.4
lambda=1.8*10^2
delta=1

x11()
D = c(1.5,2,3,5,10,100,200,400)
#D = c(5,10,15,25,50,100,200,400) #Valeurs article Baker et al, 2010

ps=0.33 #valeur Sylvain
pv=0.12 #valeur Sylvain
cv=1. #valeur Sylvain
H=0.5 #valeur Sylvain
lambda=1 #valeur Sylvain
delta=10 #valeur Sylvain

y2DM2=esp2DM2(D,ps,pv,cv,H,lambda)
se2DM2=sqrt(varmod2DM2(D,ps,pv,cv,H,lambda,delta))

plot(y2DM2~D,type="o",ylim=c(0,2))
arrows(D,y2DM2,D,y2DM2+1.96*se2DM2,length=0.05,angle=90)
arrows(D,y2DM2,D,y2DM2-1.96*se2DM2,length=0.05,angle=90)

D <- c(1.5,2,3,5,10,100,200,400)
#D <- c(1.5,2,3,5,10)

plot(D,varmod2DM2(D,ps,pv,cv,H,lambda,delta),type="o")
lines(D,varmod2DM2(D,ps,pv,cv,H,lambda,delta))

plot(D,sqrt(varmod2DM2(D,ps,pv,cv,H,lambda,delta)),type="o")


###Pour les couleurs l'indice doit être < 657

#######################################################################
#Analyse de sensibilité
D <- c(1.5,2,3,5,10,100,200,400)
#D <- c(1.5,2,3,5,10)
#D = c(5,10,15,25,50,100,200,400) #Valeurs article Baker et al, 2010

ps=0.33 #valeur Sylvain
pv=0.12 #valeur Sylvain
cv=1. #valeur Sylvain
H=5 #valeur Sylvain
lambda=1 #valeur Sylvain
delta=10 #valeur Sylvain

#-------------------------------------------------
#Sensibilité à ps
#-------------------------------------------------

ps_AS=seq(0.1,0.9,by=0.1)
y2DM2=esp2DM2(D,ps_AS[1],pv,cv,H,lambda)
se2DM2=sqrt(varmod2DM2(D,ps_AS[1],pv,cv,H,lambda,delta))

plot(y2DM2~D,type="o",pch=1,ylim=c(0,2),col=colors()[1*10],main="ps")
arrows(D,y2DM2,D,y2DM2+1.96*se2DM2,length=0.05,angle=90,col=colors()[1*10])
arrows(D,y2DM2,D,y2DM2-1.96*se2DM2,length=0.05,angle=90,col=colors()[1*10])

for (i in 2:length(ps_AS))
{
    y2DM2=esp2DM2(D,ps_AS[i],pv,cv,H,lambda)
    se2DM2=sqrt(varmod2DM2(D,ps_AS[i],pv,cv,H,lambda,delta))
    lines(y2DM2~D,type="o",pch=i)
    arrows(D,y2DM2,D,y2DM2+1.96*se2DM2,length=0.05,angle=90,col=colors()[i*10])
    arrows(D,y2DM2,D,y2DM2-1.96*se2DM2,length=0.05,angle=90,col=colors()[i*10])
}
ps=0.33

#-------------------------------------------------
#Sensibilité de E[R] et var[R] à ps; avec effets aléatoires
#-------------------------------------------------

D <- c(1.5,2,3,5,10,15,25,50,100,200,400)
#D <- c(5,10,15,25,50,100,200,400)
N <- 100 # nombre d'individus (perdrix)

#mu_ps <- log(0.33/0.67) #pour générer ps avec espérance 0.33
#sd_ps <- 0.25
#ps_p <- rnorm(N,mean=mu_ps,sd=sd_ps)
#ps <- exp(ps_p)/(1+exp(ps_p))
##ps <- rep(exp(mu_ps)/(1+exp(mu_ps)),N)
ps_nom=0.33

mu_pv <- log(0.12/0.88) #pour générer pv avec espérance 0.12
sd_pv <- 0.25
pv_p <- rnorm(N,mean=mu_pv,sd=sd_pv)
pv <- exp(pv_p)/(1+exp(pv_p))
#pv <- rep(exp(mu_pv)/(1+exp(mu_pv)),N)

mu_cv <- 1
sd_cv <- 0.25 #0.10
logmu_cv=log(mu_cv)-log(1+sd_cv^2/((mu_cv)^2))/2
logsd_cv=sqrt(log(1+sd_cv^2/((mu_cv)^2)))
cv <- exp(rnorm(N,mean=logmu_cv,sd=logsd_cv))

mu_H <- 0.5
sd_H <- 0.25
logmu_H=log(mu_H)-log(1+sd_H^2/((mu_H)^2))/2
logsd_H=sqrt(log(1+sd_H^2/((mu_H)^2)))
H <- exp(rnorm(N,mean=logmu_H,sd=logsd_H))

mu_lambda <- 1
sd_lambda <- 0.25 #0.10
logmu_lambda=log(mu_lambda)-log(1+sd_lambda^2/((mu_lambda)^2))/2
logsd_lambda=sqrt(log(1+sd_lambda^2/((mu_lambda)^2)))
lambda <- exp(rnorm(N,mean=logmu_lambda,sd=logsd_lambda))

ps_AS_seq=seq(0.1,0.9,by=0.1)
ps_AS_seq=c(ps_nom,ps_AS_seq)
mu_ps_AS <- log(ps_AS_seq/(1-ps_AS_seq))
sd_ps_AS <- 0.25

delta=10

#Génerer les val aléatoires des autres param
#Calculer moyenne et variance de la réponse fonctionnelle en même temps
x11()
y2DM2_m=matrix(0, nrow = length(mu_ps_AS), ncol = length(D))
se2DM2_m=matrix(0, nrow = length(mu_ps_AS), ncol = length(D))
for (i in 1:length(mu_ps_AS))
{
	ps_p_AS <- rnorm(N,mean=mu_ps_AS[i],sd=sd_ps_AS)
	ps_AS <- exp(ps_p_AS)/(1+exp(ps_p_AS))
	y2DM2=matrix(0, nrow = N, ncol = length(D))
	se2DM2=matrix(0, nrow = N, ncol = length(D))
for (j in 1:N)
{
    y2DM2[j,]=esp2DM2(D,ps_AS[j],pv[j],cv[j],H[j],lambda[j])
    se2DM2[j,]=sqrt(varmod2DM2(D,ps_AS[j],pv[j],cv[j],H[j],lambda[j],delta))
    }
y2DM2_m[i,] <- colMeans(y2DM2)
se2DM2_m[i,] <- colMeans(se2DM2)
#par(mfrow=c(1,2))
if (i==1) {
#plot(y2DM2_m[i,]~D,type="o",pch=1,ylim=c(0,2),col=colors()[1*9],main="ps")
plot(se2DM2_m[i,]~D,type="o",pch=1,ylim=c(0,0.3),col=colors()[1*9],main="ps")

} else {
#lines(y2DM2_m[i,]~D,type="o",pch=1,ylim=c(0,2),col=colors()[i*9],main="ps")
lines(se2DM2_m[i,]~D,type="o",pch=1,ylim=c(0,0.3),col=colors()[i*9],main="ps")
legend("topright", legend = c("Valeur nominale ps=0.33"), col = colors()[1*9], lty=1:2, cex=0.8)
}
}
ps=0.33


#-------------------------------------------------
#-------------------------------------------------
#Sensibilité à pv
#-------------------------------------------------

x11()
pv_AS=seq(0.1,0.9,by=0.1)
y2DM2=esp2DM2(D,ps,pv_AS[1],cv,H,lambda)
se2DM2=sqrt(varmod2DM2(D,ps,pv_AS[1],cv,H,lambda,delta))

plot(y2DM2~D,type="o",pch=1,ylim=c(0,1.5),col=colors()[1*10],main="pv")
arrows(D,y2DM2,D,y2DM2+1.96*se2DM2,length=0.05,angle=90,col=colors()[1*10])
arrows(D,y2DM2,D,y2DM2-1.96*se2DM2,length=0.05,angle=90,col=colors()[1*10])

for (i in 2:length(pv_AS))
{
    y2DM2=esp2DM2(D,ps,pv_AS[i],cv,H,lambda)
    se2DM2=sqrt(varmod2DM2(D,ps,pv_AS[i],cv,H,lambda,delta))
    lines(y2DM2~D,type="o",pch=i,col=colors()[i*10])
    arrows(D,y2DM2,D,y2DM2+1.96*se2DM2,length=0.05,angle=90,col=colors()[i*10])
    arrows(D,y2DM2,D,y2DM2-1.96*se2DM2,length=0.05,angle=90,col=colors()[i*10])
}
pv=0.12

#-------------------------------------------------
#Sensibilité de E[R] et var[R] à pv; avec effets aléatoires
#-------------------------------------------------

D <- c(1.5,2,3,5,10,15,25,50,100,200,400)
N <- 100 # nombre d'individus (perdrix)

mu_ps <- log(0.33/0.67) #pour générer ps avec espérance 0.33
sd_ps <- 0.25
ps_p <- rnorm(N,mean=mu_ps,sd=sd_ps)
ps <- exp(ps_p)/(1+exp(ps_p))
#ps <- rep(exp(mu_ps)/(1+exp(mu_ps)),N)

#mu_pv <- log(0.12/0.88) #pour générer pv avec espérance 0.12
#sd_pv <- 0.25
#pv_p <- rnorm(N,mean=mu_pv,sd=sd_pv)
#pv <- exp(pv_p)/(1+exp(pv_p))
##pv <- rep(exp(mu_pv)/(1+exp(mu_pv)),N)
pv_nom=0.12

mu_cv <- 1
sd_cv <- 0.25 #0.10
logmu_cv=log(mu_cv)-log(1+sd_cv^2/((mu_cv)^2))/2
logsd_cv=sqrt(log(1+sd_cv^2/((mu_cv)^2)))
cv <- exp(rnorm(N,mean=logmu_cv,sd=logsd_cv))

mu_H <- 0.5
sd_H <- 0.25
logmu_H=log(mu_H)-log(1+sd_H^2/((mu_H)^2))/2
logsd_H=sqrt(log(1+sd_H^2/((mu_H)^2)))
H <- exp(rnorm(N,mean=logmu_H,sd=logsd_H))

mu_lambda <- 1
sd_lambda <- 0.25 #0.10
logmu_lambda=log(mu_lambda)-log(1+sd_lambda^2/((mu_lambda)^2))/2
logsd_lambda=sqrt(log(1+sd_lambda^2/((mu_lambda)^2)))
lambda <- exp(rnorm(N,mean=logmu_lambda,sd=logsd_lambda))


pv_AS_seq=seq(0.1,0.9,by=0.1)
pv_AS_seq=c(pv_nom,pv_AS_seq)
mu_pv_AS <- log(pv_AS_seq/(1-pv_AS_seq))
sd_pv_AS <- 0.25

delta=10

#Génerer les val aléatoires des autres param
#Calculer moyenne et variance de la réponse fonctionnelle en même temps
x11()
y2DM2_m=matrix(0, nrow = length(mu_pv_AS), ncol = length(D))
se2DM2_m=matrix(0, nrow = length(mu_pv_AS), ncol = length(D))
for (i in 1:length(mu_pv_AS))
{
	pv_p_AS <- rnorm(N,mean=mu_pv_AS[i],sd=sd_pv_AS)
	pv_AS <- exp(pv_p_AS)/(1+exp(pv_p_AS))
	
	y2DM2=matrix(0, nrow = N, ncol = length(D))
	se2DM2=matrix(0, nrow = N, ncol = length(D))
for (j in 1:N)
{
    y2DM2[j,]=esp2DM2(D,ps[j],pv_AS[j],cv[j],H[j],lambda[j])
    se2DM2[j,]=sqrt(varmod2DM2(D,ps[j],pv_AS[j],cv[j],H[j],lambda[j],delta))
    }
y2DM2_m[i,] <- colMeans(y2DM2)
se2DM2_m[i,] <- colMeans(se2DM2)
#par(mfrow=c(1,2))
if (i==1) {
plot(y2DM2_m[i,]~D,type="o",pch=1,ylim=c(0,2),col=colors()[1*9],main="pv")
#plot(se2DM2_m[i,]~D,type="o",pch=1,ylim=c(0,0.3),col=colors()[1*9],main="pv")
legend("topright", legend = c("Valeur nominale pv=0.12"), col = colors()[1*9], lty=1:2, cex=0.8)

} else {
lines(y2DM2_m[i,]~D,type="o",pch=1,ylim=c(0,2),col=colors()[i*9],main="pv")
#lines(se2DM2_m[i,]~D,type="o",pch=1,ylim=c(0,0.3),col=colors()[i*9],main="pv")
}
}
pv=0.12

#_________________________________
#_________________________________
#Sensibilité à cv
#-------------------------------------------------

x11()
D <- c(1.5,2,3,5,10,15,25,50,100,200,400)
y2DM2=esp2DM2(D,ps,pv,cv_AS[1],H,lambda)
se2DM2=sqrt(varmod2DM2(D,ps,pv,cv_AS[1],H,lambda,delta))

plot(y2DM2~D,type="o",pch=1,ylim=c(0,1.5),col=colors()[1*10],main="cv")
arrows(D,y2DM2,D,y2DM2+1.96*se2DM2,length=0.05,angle=90,col=colors()[1*10])
arrows(D,y2DM2,D,y2DM2-1.96*se2DM2,length=0.05,angle=90,col=colors()[1*10])

for (i in 2:length(cv_AS))
{
    y2DM2=esp2DM2(D,ps,pv,cv_AS[i],H,lambda)
    se2DM2=sqrt(varmod2DM2(D,ps,pv,cv_AS[i],H,lambda,delta))
    lines(y2DM2~D,type="o",pch=i,col=colors()[i*10])
    arrows(D,y2DM2,D,y2DM2+1.96*se2DM2,length=0.05,angle=90,col=colors()[i*10])
    arrows(D,y2DM2,D,y2DM2-1.96*se2DM2,length=0.05,angle=90,col=colors()[i*10])
}
cv=1

#-----------------------------------------
#Sensibilité de E[R] et var[R] à cv; avec effets aléatoires
#-----------------------------------------


D <- c(1.5,2,3,5,10,15,25,50,100,200,400)
N <- 100 # nombre d'individus (perdrix)

mu_ps <- log(0.33/0.67) #pour générer ps avec espérance 0.33
sd_ps <- 0.25
ps_p <- rnorm(N,mean=mu_ps,sd=sd_ps)
ps <- exp(ps_p)/(1+exp(ps_p))
##ps <- rep(exp(mu_ps)/(1+exp(mu_ps)),N)

mu_pv <- log(0.12/0.88) #pour générer pv avec espérance 0.12
sd_pv <- 0.25
pv_p <- rnorm(N,mean=mu_pv,sd=sd_pv)
pv <- exp(pv_p)/(1+exp(pv_p))
#pv <- rep(exp(mu_pv)/(1+exp(mu_pv)),N)

#mu_cv <- 1
#sd_cv <- 0.25 #0.10
#logmu_cv=log(mu_cv)-log(1+sd_cv^2/((mu_cv)^2))/2
#logsd_cv=sqrt(log(1+sd_cv^2/((mu_cv)^2)))
#cv <- exp(rnorm(N,mean=logmu_cv,sd=logsd_cv))
cv_nom=1

mu_H <- 0.5
sd_H <- 0.25
logmu_H=log(mu_H)-log(1+sd_H^2/((mu_H)^2))/2
logsd_H=sqrt(log(1+sd_H^2/((mu_H)^2)))
H <- exp(rnorm(N,mean=logmu_H,sd=logsd_H))

mu_lambda <- 1
sd_lambda <- 0.25 #0.10
logmu_lambda=log(mu_lambda)-log(1+sd_lambda^2/((mu_lambda)^2))/2
logsd_lambda=sqrt(log(1+sd_lambda^2/((mu_lambda)^2)))
lambda <- exp(rnorm(N,mean=logmu_lambda,sd=logsd_lambda))


mu_cv_AS=c(0.3,2,5,10)
mu_cv_AS=c(cv_nom,mu_cv_AS)
sd_cv_AS <- 0.25
logmu_cv_AS=log(mu_cv_AS)-log(1+sd_cv_AS^2/((mu_cv_AS)^2))/2
logsd_cv_AS=sqrt(log(1+sd_cv_AS^2/((mu_cv_AS)^2)))

delta=10

#Génerer les val aléatoires des autres param
#Calculer moyenne et variance de la réponse fonctionnelle en même temps
x11()
se2DM2_m=matrix(0, nrow = length(mu_cv_AS), ncol = length(D))
for (i in 1:length(mu_cv_AS))
{
	#Attention appliquer la bonne transformaton entre la normale et la log-normale	
	cv_AS <- exp(rnorm(N,mean=logmu_cv_AS[i],sd=logsd_cv_AS))
	se2DM2=matrix(0, nrow = N, ncol = length(D))
for (j in 1:N)
{   
    se2DM2[j,]=sqrt(varmod2DM2(D,ps[j],pv[j],cv_AS[j],H[j],lambda[j],delta))
    }
se2DM2_m[i,] <- colMeans(se2DM2)

if (i==1) {
plot(se2DM2_m[i,]~D,type="o",pch=1,ylim=c(0,0.5),col=colors()[i*9],main="cv")
legend("topright", legend = c("Valeur nominale cv=1"), col = colors()[1*9], lty=1:2, cex=0.8)

} else {
lines(se2DM2_m[i,]~D,type="o",pch=1,ylim=c(0,0.1),col=colors()[i*9],main="cv")
}
}
cv=1

#-----------------------------------------
#-----------------------------------------
#Sensibilité à H
#-----------------------------------------

x11()
H_AS=c(0.5,2,5,7,10)
y2DM2=esp2DM2(D,ps,pv,cv,H_AS[1],lambda)
se2DM2=sqrt(varmod2DM2(D,ps,pv,cv,H_AS[1],lambda,delta))

plot(y2DM2~D,type="o",pch=1,ylim=c(0,1.5),col=colors()[1*10],main="H")
arrows(D,y2DM2,D,y2DM2+1.96*se2DM2,length=0.05,angle=90,col=colors()[1*10])
arrows(D,y2DM2,D,y2DM2-1.96*se2DM2,length=0.05,angle=90,col=colors()[1*10])

for (i in 2:length(H_AS))
{
    y2DM2=esp2DM2(D,ps,pv,cv,H_AS[i],lambda)
    se2DM2=sqrt(varmod2DM2(D,ps,pv,cv,H_AS[i],lambda,delta))
    lines(y2DM2~D,type="o",pch=i,col=colors()[i*10])
    arrows(D,y2DM2,D,y2DM2+1.96*se2DM2,length=0.05,angle=90,col=colors()[i*10])
    arrows(D,y2DM2,D,y2DM2-1.96*se2DM2,length=0.05,angle=90,col=colors()[i*10])
}
H=0.5

#-----------------------------------------
#Sensibilité de E[R] et var[R] à H; avec effets aléatoires
#-----------------------------------------

D <- c(1.5,2,3,5,10,15,25,50,100,200,400)
N <- 100 # nombre d'individus (perdrix)

mu_ps <- log(0.33/0.67) #pour générer ps avec espérance 0.33
sd_ps <- 0.25
ps_p <- rnorm(N,mean=mu_ps,sd=sd_ps)
ps <- exp(ps_p)/(1+exp(ps_p))
##ps <- rep(exp(mu_ps)/(1+exp(mu_ps)),N)

mu_pv <- log(0.12/0.88) #pour générer pv avec espérance 0.12
sd_pv <- 0.25
pv_p <- rnorm(N,mean=mu_pv,sd=sd_pv)
pv <- exp(pv_p)/(1+exp(pv_p))
#pv <- rep(exp(mu_pv)/(1+exp(mu_pv)),N)

mu_cv <- 1
sd_cv <- 0.25 #0.01
logmu_cv=log(mu_cv)-log(1+sd_cv^2/((mu_cv)^2))/2
logsd_cv=sqrt(log(1+sd_cv^2/((mu_cv)^2)))
cv <- exp(rnorm(N,mean=logmu_cv,sd=logsd_cv))

#mu_H <- 0.5
#sd_H <- 0.25
#logmu_H=log(mu_H)-log(1+sd_H^2/((mu_H)^2))/2
#logsd_H=sqrt(log(1+sd_H^2/((mu_H)^2)))
#H <- exp(rnorm(N,mean=logmu_H,sd=logsd_H))
H_nom=0.5

mu_lambda <- 1
sd_lambda <- 0.25 #0.10
logmu_lambda=log(mu_lambda)-log(1+sd_lambda^2/((mu_lambda)^2))/2
logsd_lambda=sqrt(log(1+sd_lambda^2/((mu_lambda)^2)))
lambda <- exp(rnorm(N,mean=logmu_lambda,sd=logsd_lambda))


mu_H_AS=c(2,5,7,10)
mu_H_AS=c(H_nom,mu_H_AS)
sd_H_AS <- 0.25
logmu_H_AS=log(mu_H_AS)-log(1+sd_H_AS^2/((mu_H_AS)^2))/2
logsd_H_AS=sqrt(log(1+sd_H_AS^2/((mu_H_AS)^2)))

delta=10

#Génerer les val aléatoires des autres param
#Calculer moyenne et variance de la réponse fonctionnelle en même temps
x11()
y2DM2_m=matrix(0, nrow = length(mu_H_AS), ncol = length(D))
se2DM2_m=matrix(0, nrow = length(mu_H_AS), ncol = length(D))
for (i in 1:length(mu_H_AS))
{
	#Attention appliquer la bonne transformaton entre la normale et la log-normale	
	H_AS <- exp(rnorm(N,mean=logmu_H_AS[i],sd=logsd_H_AS))
	
	y2DM2=matrix(0, nrow = N, ncol = length(D))
	se2DM2=matrix(0, nrow = N, ncol = length(D))
for (j in 1:N)
{   
	y2DM2[j,]=esp2DM2(D,ps[j],pv[j],cv[j],H_AS[j],lambda[j])
    se2DM2[j,]=sqrt(varmod2DM2(D,ps[j],pv[j],cv[j],H_AS[j],lambda[j],delta))
    }
y2DM2_m[i,] <- colMeans(y2DM2)
se2DM2_m[i,] <- colMeans(se2DM2)

if (i==1) {
plot(y2DM2_m[i,]~D,type="o",pch=1,ylim=c(0,1.5),col=colors()[1*9],main="H")
#plot(se2DM2_m[i,]~D,type="o",pch=1,ylim=c(0,0.3),col=colors()[i*9],main="H")
legend("topright", legend = c("Valeur nominale H=5"), col = colors()[1*9], lty=1:2, cex=0.8)

} else {
lines(y2DM2_m[i,]~D,type="o",pch=1,ylim=c(0,1.5),col=colors()[i*9],main="H")
#lines(se2DM2_m[i,]~D,type="o",pch=1,ylim=c(0,0.3),col=colors()[i*9],main="H")
}
}
H=0.5

#-----------------------------------------
#-----------------------------------------
#Sensibilité à lambda
#-----------------------------------------

x11()
lambda_AS=c(1,2,5,10,50,100,200)
y2DM2=esp2DM2(D,ps,pv,cv,H,lambda_AS[1])
se2DM2=sqrt(varmod2DM2(D,ps,pv,cv,H,lambda_AS[1],delta))

plot(y2DM2~D,type="o",pch=1,ylim=c(0,1.5),col=colors()[1*10],main="lambda")
arrows(D,y2DM2,D,y2DM2+1.96*se2DM2,length=0.05,angle=90,col=colors()[1*10])
arrows(D,y2DM2,D,y2DM2-1.96*se2DM2,length=0.05,angle=90,col=colors()[1*10])

for (i in 2:length(lambda_AS))
{
    y2DM2=esp2DM2(D,ps,pv,cv,H,lambda_AS[i])
    se2DM2=sqrt(varmod2DM2(D,ps,pv,cv,H,lambda_AS[i],delta))
    lines(y2DM2~D,type="o",pch=i,col=colors()[i*10])
    arrows(D,y2DM2,D,y2DM2+1.96*se2DM2,length=0.05,angle=90,col=colors()[i*10])
    arrows(D,y2DM2,D,y2DM2-1.96*se2DM2,length=0.05,angle=90,col=colors()[i*10])
}
lambda=1

#-----------------------------------------
#Sensibilité de E[R] et var[R] à lambda; avec effets aléatoires
#-----------------------------------------

D <- c(1.5,2,3,5,10,15,25,50,100,200,400)
N <- 100 # nombre d'individus (perdrix)

mu_ps <- log(0.33/0.67) #pour générer ps avec espérance 0.33
sd_ps <- 0.25
ps_p <- rnorm(N,mean=mu_ps,sd=sd_ps)
ps <- exp(ps_p)/(1+exp(ps_p))
##ps <- rep(exp(mu_ps)/(1+exp(mu_ps)),N)

mu_pv <- log(0.12/0.88) #pour générer pv avec espérance 0.12
sd_pv <- 0.25
pv_p <- rnorm(N,mean=mu_pv,sd=sd_pv)
pv <- exp(pv_p)/(1+exp(pv_p))
#pv <- rep(exp(mu_pv)/(1+exp(mu_pv)),N)

mu_cv <- 1
sd_cv <- 0.25 #0.01
logmu_cv=log(mu_cv)-log(1+sd_cv^2/((mu_cv)^2))/2
logsd_cv=sqrt(log(1+sd_cv^2/((mu_cv)^2)))
cv <- exp(rnorm(N,mean=logmu_cv,sd=logsd_cv))

mu_H <- 0.5
sd_H <- 0.25
logmu_H=log(mu_H)-log(1+sd_H^2/((mu_H)^2))/2
logsd_H=sqrt(log(1+sd_H^2/((mu_H)^2)))
H <- exp(rnorm(N,mean=logmu_H,sd=logsd_H))

#mu_lambda <- 1
#sd_lambda <- 0.1 #0.10
#logmu_lambda=log(mu_lambda)-log(1+sd_lambda^2/((mu_lambda)^2))/2
#logsd_lambda=sqrt(log(1+sd_lambda^2/((mu_lambda)^2)))
#lambda <- exp(rnorm(N,mean=logmu_lambda,sd=logsd_lambda))
lambda_num=1


mu_lambda_AS=c(0.1,2,5,10,100)
mu_lambda_AS=c(lambda_num,mu_lambda_AS)
sd_lambda_AS <- 0.25
logmu_lambda_AS=log(mu_lambda_AS)-log(1+sd_lambda_AS^2/((mu_lambda_AS)^2))/2
logsd_lambda_AS=sqrt(log(1+sd_lambda_AS^2/((mu_lambda_AS)^2)))

delta=10

#Génerer les val aléatoires des autres param
#Calculer moyenne et variance de la réponse fonctionnelle en même temps
x11()
y2DM2_m=matrix(0, nrow = length(mu_lambda_AS), ncol = length(D))
se2DM2_m=matrix(0, nrow = length(mu_lambda_AS), ncol = length(D))
for (i in 1:length(mu_lambda_AS))
{
	#Attention appliquer la bonne transformaton entre la normale et la log-normale	
	lambda_AS <- exp(rnorm(N,mean=logmu_lambda_AS[i],sd=logsd_lambda_AS))
	
	y2DM2=matrix(0, nrow = N, ncol = length(D))
	se2DM2=matrix(0, nrow = N, ncol = length(D))
for (j in 1:N)
{   
	y2DM2[j,]=esp2DM2(D,ps[j],pv[j],cv[j],H[j],lambda_AS[j])
    se2DM2[j,]=sqrt(varmod2DM2(D,ps[j],pv[j],cv[j],H[j],lambda_AS[j],delta))
    }
y2DM2_m[i,] <- colMeans(y2DM2)
se2DM2_m[i,] <- colMeans(se2DM2)

if (i==1) {
plot(y2DM2_m[i,]~D,type="o",pch=1,ylim=c(0,2),col=colors()[1*9],main="lambda")
#plot(se2DM2_m[i,]~D,type="o",pch=1,ylim=c(0,0.5),col=colors()[i*9],main="lambda")
legend("topright", legend = c("Valeur nominale lambda=1"), col = colors()[1*9], lty=1:2, cex=0.8)

} else {
lines(y2DM2_m[i,]~D,type="o",pch=1,ylim=c(0,2),col=colors()[i*9],main="lambda")
#lines(se2DM2_m[i,]~D,type="o",pch=1,ylim=c(0,0.5),col=colors()[i*9],main="lambda")
}
}
lambda=1

#-----------------------------------------
#-----------------------------------------
#Sensibilité à Delta
#-------------------------------------------------

x11()
delta_AS=c(1,10,100,1000)
y2DM2=esp2DM2(D,ps,pv,cv,H,lambda)
se2DM2=sqrt(varmod2DM2(D,ps,pv,cv,H,lambda,delta_AS[1]))

plot(y2DM2~D,type="o",pch=1,ylim=c(0,2),col=colors()[1*10],main="Delta")
arrows(D,y2DM2,D,y2DM2+1.96*se2DM2,length=0.05,angle=90,col=colors()[1*10])
arrows(D,y2DM2,D,y2DM2-1.96*se2DM2,length=0.05,angle=90,col=colors()[1*10])

for (i in 2:length(delta_AS))
{
    y2DM2=esp2DM2(D,ps,pv,cv,H,lambda)
    se2DM2=sqrt(varmod2DM2(D,ps,pv,cv,H,lambda,delta_AS[i]))
    lines(y2DM2~D,type="o",pch=i,col=colors()[i*10])
    arrows(D,y2DM2,D,y2DM2+1.96*se2DM2,length=0.05,angle=90,col=colors()[i*10])
    arrows(D,y2DM2,D,y2DM2-1.96*se2DM2,length=0.05,angle=90,col=colors()[i*10])
}
delta=10

#-----------------------------------------
#Sensibilité de E[R] et var[R] à delta; avec effets aléatoires
#-----------------------------------------

D <- c(1.5,2,3,5,10,15,25,50,100,200,400)
N <- 100 # nombre d'individus (perdrix)

mu_ps <- log(0.33/0.67) #pour générer ps avec espérance 0.33
sd_ps <- 0.25
ps_p <- rnorm(N,mean=mu_ps,sd=sd_ps)
ps <- exp(ps_p)/(1+exp(ps_p))
##ps <- rep(exp(mu_ps)/(1+exp(mu_ps)),N)

mu_pv <- log(0.12/0.88) #pour générer pv avec espérance 0.12
sd_pv <- 0.25
pv_p <- rnorm(N,mean=mu_pv,sd=sd_pv)
pv <- exp(pv_p)/(1+exp(pv_p))
#pv <- rep(exp(mu_pv)/(1+exp(mu_pv)),N)

mu_cv <- 1
sd_cv <- 0.25 #0.01
logmu_cv=log(mu_cv)-log(1+sd_cv^2/((mu_cv)^2))/2
logsd_cv=sqrt(log(1+sd_cv^2/((mu_cv)^2)))
cv <- exp(rnorm(N,mean=logmu_cv,sd=logsd_cv))

mu_H <- 0.5
sd_H <- 0.25
logmu_H=log(mu_H)-log(1+sd_H^2/((mu_H)^2))/2
logsd_H=sqrt(log(1+sd_H^2/((mu_H)^2)))
H <- exp(rnorm(N,mean=logmu_H,sd=logsd_H))

mu_lambda <- 1
sd_lambda <- 0.25 #0.10
logmu_lambda=log(mu_lambda)-log(1+sd_lambda^2/((mu_lambda)^2))/2
logsd_lambda=sqrt(log(1+sd_lambda^2/((mu_lambda)^2)))
lambda <- exp(rnorm(N,mean=logmu_lambda,sd=logsd_lambda))

delta_num=10

delta_AS=c(1,5,100)
delta_AS=c(delta_num,delta_AS)

#Génerer les val aléatoires des autres param
#Calculer moyenne et variance de la réponse fonctionnelle en même temps
#y2DM2_m=matrix(0, nrow = length(mu_delta_AS), ncol = length(D))
x11()
se2DM2_m=matrix(0, nrow = length(delta_AS), ncol = length(D))
for (i in 1:length(delta_AS))
{
	#Attention appliquer la bonne transformaton entre la normale et la log-normale	
	
#	y2DM2=matrix(0, nrow = N, ncol = length(D))
	se2DM2=matrix(0, nrow = N, ncol = length(D))
for (j in 1:N)
{   
#	y2DM2[j,]=esp2DM2(D,ps[j],pv[j],cv[j],H[j],lambda[j])
    se2DM2[j,]=sqrt(varmod2DM2(D,ps[j],pv[j],cv[j],H[j],lambda[j],delta_AS[i]))
    }
#y2DM2_m[i,] <- colMeans(y2DM2)
se2DM2_m[i,] <- colMeans(se2DM2)

if (i==1) {
plot(se2DM2_m[i,]~D,type="o",pch=1,ylim=c(0,0.2),col=colors()[i*9],main="delta")
legend("topright", legend = c("Valeur nominale delta=10"), col = colors()[1*9], lty=1:2, cex=0.8)

} else {
lines(se2DM2_m[i,]~D,type="o",pch=1,ylim=c(0,0.2),col=colors()[i*9],main="delta")
}
}
delta=10

#-----------------------------------------
#-----------------------------------------




esp_l=lambda; var_l=0.05
esp_logl=log(esp_l)-0.5*log(1+var_l/((esp_l)^2))
var_logl=log(1+var_l/((esp_l)^2))
sd_logl=sqrt(var_logl)

esp_H=H; var_H=0.1
esp_logH=log(esp_H)-0.5*log(1+var_H/((esp_H)^2))
var_logH=log(1+var_H/((esp_H)^2))
sd_logH=sqrt(var_logH)

for (i in 2:5)
{
    log_l=rnorm(1,esp_logl,sd_logl)
    lambda=exp(log_l)
    log_H=rnorm(1,esp_logH,sd_logH)
    H=exp(log_H)
    par(new=T)
    y2D=esp2D(lambda,D,H)
    se2D=sqrt(varmod2D(lambda,D,H,delta))
    lines(y2D~D,type="o",pch=i)
    arrows(D,y2D,D,y2D+1.96*se2D,length=0.05,angle=90)
    arrows(D,y2D,D,y2D-1.96*se2D,length=0.05,angle=90)
}




##Couleurs
plot(1:27, type = "n", axes = FALSE,
xlab = "", ylab = "", main = "Couleurs sous R (col=)")
k <- 0 # k, indice du vecteur des couleurs à parcourir
for(i in 1:26) { # on fait varier les lignes
    for(j in 1:26) { # on fait varier les colonnes
        k <- k + 1
        if (k > 657)
        break # Taille du vecteur à ne pas dépasser
        polygon(c(i, i+1, i+1, i, i),c(j, j, j+1, j+1, j),
        col = colors()[k])
        text(i+0.5, j+0.5, as.character(k), cex = 0.8)
    }
}

#legend("topright", 
#  legend = c("Group 1", "Group 2"), 
#  col = c(rgb(0.2,0.4,0.1,0.7), 
#  rgb(0.8,0.4,0.1,0.7)), 
#  pch = c(17,19), 
#  bty = "n", 
#  pt.cex = 2, 
#  cex = 1.2, 
#  text.col = "black", 
#  horiz = F , 
#  inset = c(0.1, 0.1))

}
#-----------------------------------------
#-----------------------------------------

#Représentation de secénarios contrastés

D <- c(1.5,2,3,5,10,15,25,50,100,200,400)
N <- 100 # nombre d'individus (perdrix)

#Scénario dit de Sylvain
mu_ps <- log(0.33/0.67) #pour générer ps avec espérance 0.33
sd_ps <- 0.25
ps_p <- rnorm(N,mean=mu_ps,sd=sd_ps)
ps <- exp(ps_p)/(1+exp(ps_p))

mu_pv <- log(0.12/0.88) #pour générer pv avec espérance 0.12
sd_pv <- 0.25
pv_p <- rnorm(N,mean=mu_pv,sd=sd_pv)
pv <- exp(pv_p)/(1+exp(pv_p))

mu_cv <- 1
sd_cv <- 0.25 #0.10
logmu_cv=log(mu_cv)-log(1+sd_cv^2/((mu_cv)^2))/2
logsd_cv=sqrt(log(1+sd_cv^2/((mu_cv)^2)))
cv <- exp(rnorm(N,mean=logmu_cv,sd=logsd_cv))

mu_H <- 0.5
sd_H <- 0.25
logmu_H=log(mu_H)-log(1+sd_H^2/((mu_H)^2))/2
logsd_H=sqrt(log(1+sd_H^2/((mu_H)^2)))
H <- exp(rnorm(N,mean=logmu_H,sd=logsd_H))

mu_lambda <- 1
sd_lambda <- 0.25 #0.10
logmu_lambda=log(mu_lambda)-log(1+sd_lambda^2/((mu_lambda)^2))/2
logsd_lambda=sqrt(log(1+sd_lambda^2/((mu_lambda)^2)))
lambda <- exp(rnorm(N,mean=logmu_lambda,sd=logsd_lambda))

delta=10

#Scénario avec espérance et variance de réponse fonctionnelle plus faibles
mu_ps <- log(0.33/0.67)
sd_ps <- 0.25
ps_p <- rnorm(N,mean=mu_ps,sd=sd_ps)
ps <- exp(ps_p)/(1+exp(ps_p))

mu_pv <- log(0.12/0.88)
sd_pv <- 0.25
pv_p <- rnorm(N,mean=mu_pv,sd=sd_pv)
pv <- exp(pv_p)/(1+exp(pv_p))

mu_cv <- 1
sd_cv <- 0.25 #0.10
logmu_cv=log(mu_cv)-log(1+sd_cv^2/((mu_cv)^2))/2
logsd_cv=sqrt(log(1+sd_cv^2/((mu_cv)^2)))
cv <- exp(rnorm(N,mean=logmu_cv,sd=logsd_cv))

mu_H <- 5
sd_H <- 0.25
logmu_H=log(mu_H)-log(1+sd_H^2/((mu_H)^2))/2
logsd_H=sqrt(log(1+sd_H^2/((mu_H)^2)))
H <- exp(rnorm(N,mean=logmu_H,sd=logsd_H))

mu_lambda <- 1
sd_lambda <- 0.25 #0.10
logmu_lambda=log(mu_lambda)-log(1+sd_lambda^2/((mu_lambda)^2))/2
logsd_lambda=sqrt(log(1+sd_lambda^2/((mu_lambda)^2)))
lambda <- exp(rnorm(N,mean=logmu_lambda,sd=logsd_lambda))

delta=10

#Scénario avec espérance de réponse fonctionnelle un peu plus basse et variance aussi
mu_ps <- log(0.9/0.1)
sd_ps <- 0.25
ps_p <- rnorm(N,mean=mu_ps,sd=sd_ps)
ps <- exp(ps_p)/(1+exp(ps_p))

mu_pv <- log(0.12/0.88)
sd_pv <- 0.25
pv_p <- rnorm(N,mean=mu_pv,sd=sd_pv)
pv <- exp(pv_p)/(1+exp(pv_p))

mu_cv <- 1
sd_cv <- 0.25 #0.10
logmu_cv=log(mu_cv)-log(1+sd_cv^2/((mu_cv)^2))/2
logsd_cv=sqrt(log(1+sd_cv^2/((mu_cv)^2)))
cv <- exp(rnorm(N,mean=logmu_cv,sd=logsd_cv))

mu_H <- 1
sd_H <- 0.25
logmu_H=log(mu_H)-log(1+sd_H^2/((mu_H)^2))/2
logsd_H=sqrt(log(1+sd_H^2/((mu_H)^2)))
H <- exp(rnorm(N,mean=logmu_H,sd=logsd_H))

mu_lambda <- 1
sd_lambda <- 0.25 #0.10
logmu_lambda=log(mu_lambda)-log(1+sd_lambda^2/((mu_lambda)^2))/2
logsd_lambda=sqrt(log(1+sd_lambda^2/((mu_lambda)^2)))
lambda <- exp(rnorm(N,mean=logmu_lambda,sd=logsd_lambda))

delta=10


#Scénario avec espérance de réponse fonctionnelle un peu plus basse et variance aussi
mu_ps <- log(0.9/0.1)
sd_ps <- 0.25
ps_p <- rnorm(N,mean=mu_ps,sd=sd_ps)
ps <- exp(ps_p)/(1+exp(ps_p))

mu_pv <- log(0.12/0.88)
sd_pv <- 0.25
pv_p <- rnorm(N,mean=mu_pv,sd=sd_pv)
pv <- exp(pv_p)/(1+exp(pv_p))

mu_cv <- 1
sd_cv <- 0.25 #0.10
logmu_cv=log(mu_cv)-log(1+sd_cv^2/((mu_cv)^2))/2
logsd_cv=sqrt(log(1+sd_cv^2/((mu_cv)^2)))
cv <- exp(rnorm(N,mean=logmu_cv,sd=logsd_cv))

mu_H <- 1
sd_H <- 0.25
logmu_H=log(mu_H)-log(1+sd_H^2/((mu_H)^2))/2
logsd_H=sqrt(log(1+sd_H^2/((mu_H)^2)))
H <- exp(rnorm(N,mean=logmu_H,sd=logsd_H))

mu_lambda <- 1
sd_lambda <- 0.25 #0.10
logmu_lambda=log(mu_lambda)-log(1+sd_lambda^2/((mu_lambda)^2))/2
logsd_lambda=sqrt(log(1+sd_lambda^2/((mu_lambda)^2)))
lambda <- exp(rnorm(N,mean=logmu_lambda,sd=logsd_lambda))

delta=10


#Scénario avec variance de réponse fonctionnelle plus basse
mu_ps <- log(0.9/0.1)
sd_ps <- 0.25
ps_p <- rnorm(N,mean=mu_ps,sd=sd_ps)
ps <- exp(ps_p)/(1+exp(ps_p))

mu_pv <- log(0.12/0.88)
sd_pv <- 0.25
pv_p <- rnorm(N,mean=mu_pv,sd=sd_pv)
pv <- exp(pv_p)/(1+exp(pv_p))

mu_cv <- 1
sd_cv <- 0.25 #0.10
logmu_cv=log(mu_cv)-log(1+sd_cv^2/((mu_cv)^2))/2
logsd_cv=sqrt(log(1+sd_cv^2/((mu_cv)^2)))
cv <- exp(rnorm(N,mean=logmu_cv,sd=logsd_cv))

mu_H <- 1
sd_H <- 0.25
logmu_H=log(mu_H)-log(1+sd_H^2/((mu_H)^2))/2
logsd_H=sqrt(log(1+sd_H^2/((mu_H)^2)))
H <- exp(rnorm(N,mean=logmu_H,sd=logsd_H))

mu_lambda <- 1
sd_lambda <- 0.25 #0.10
logmu_lambda=log(mu_lambda)-log(1+sd_lambda^2/((mu_lambda)^2))/2
logsd_lambda=sqrt(log(1+sd_lambda^2/((mu_lambda)^2)))
lambda <- exp(rnorm(N,mean=logmu_lambda,sd=logsd_lambda))

delta=50


#Calculer moyenne et variance de la réponse fonctionnelle en même temps
y2DM2=matrix(0, nrow = N, ncol = length(D))
se2DM2=matrix(0, nrow = N, ncol = length(D))
for (j in 1:N)
{
    y2DM2[j,]=esp2DM2(D,ps[j],pv[j],cv[j],H[j],lambda[j])
    se2DM2[j,]=sqrt(varmod2DM2(D,ps[j],pv[j],cv[j],H[j],lambda[j],delta))
}
y2DM2_m <- colMeans(y2DM2)
se2DM2_m <- colMeans(se2DM2)

x11()
plot(y2DM2_m~D,type="o",ylim=c(0,2),col="red")
#lines(y2DM2_m~D,type="o",ylim=c(0,2),col="red")
arrows(D,y2DM2_m,D,y2DM2_m+1.96*se2DM2_m,length=0.05,angle=90)
arrows(D,y2DM2_m,D,y2DM2_m-1.96*se2DM2_m,length=0.05,angle=90)


