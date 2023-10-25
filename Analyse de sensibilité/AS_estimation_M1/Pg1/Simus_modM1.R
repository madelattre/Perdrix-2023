## Maud Delattre
## Simulations de données de réponses fonctionnelles selon le modèle (M1)

## EspR : espérance de la réponse fonctionnelle
## VarR : variance de la réponse fonctionnelle

## Plusieurs paramètres : 
## - C2E et C2V sont des constantes connues
## - lambda (L/2v)
## - ch : temps de manipulation
## - y  : densité de graines 
## - delta

EspR <- function(lambda, ch, y){
  C2E <-  (sqrt(2)+log(1+sqrt(2)))/3
  res <- (sqrt(y)-1)/(C2E*lambda+ch*(sqrt(y)-1))
  return(res)
}

VarR <- function(lambda, ch, y, delta){
  C2E <- (sqrt(2)+log(1+sqrt(2)))/3
  C2V <- 2/3
  res <- 1/delta*(C2V-(C2E)^2)*lambda^2*(sqrt(y)-1)/(C2E*lambda+ch*(sqrt(y)-1))^3
  return(res)
}