library(multisensi)

#On fait varier l'espérence et le coefficient de variation (std/esp)

ASdesign_M1=expand.grid(esp_lambda=c(1,5,10,20), var_lambda=c(0.25,0.5,0.75,1), esp_ch=c(0.5,2,5,10), var_ch=c(0.25,0.5,0.75,1))

ASdesign_M1$var_lambda=(ASdesign_M1$esp_lambda*ASdesign_M1$var_lambda)^2
ASdesign_M1$var_ch=(ASdesign_M1$esp_ch*ASdesign_M1$var_ch)^2

write.table(ASdesign_M1,'ASdesign_M1.csv',row.names=F,sep="\t")
