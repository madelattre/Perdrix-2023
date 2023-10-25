################################################################
## Modèle en 2D

CE <- (sqrt(2)+log(1+sqrt(2)))/3
CV <- 2/3

esp2DM2<-function(D,ps,pv,cv,H,lambda){
  return((ps*(1-pv)*(sqrt(D)-1))/(pv*CV*(sqrt(D)-1)+CE*lambda*(1-pv)+H*ps*(1-pv)*(sqrt(D)-1)))
}


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
H=5 #valeur Sylvain
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

#Sensibilité à ps
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


#Sensibilité à pv
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


#Sensibilité à cv
x11()
cv_AS=c(0.3,1,2,5,10)
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


#Sensibilité à H
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


#Sensibilité à lambda
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


#Sensibilité à Delta
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
