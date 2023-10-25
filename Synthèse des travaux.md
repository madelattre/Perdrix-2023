## Modèles


$Y_{ij}$ : réponse fonctionnelle de la $i$e perdrix pour une densité de graines $d_{ij}$
$N$: nombre d'individus


1. **Modèle (M1)**
$$

Y_{ij} | \psi_i = (\lambda_i,h_i)\sim \mathcal{N}(m_{ij},s_{ij}^2)

$$ 

où $$m_{ij} = \frac{\sqrt{d_{ij}}-1}{C_{2E} \textcolor{blue}{\lambda_i} + \textcolor{blue}{h_i} (\sqrt{d_{ij}}-1)}$$ et $$s_{ij}^2 = \frac{C_{2V}-C_{2E}^2}{\Delta} \frac{\textcolor{blue}{\lambda_i^2} (\sqrt{d_{ij}}-1)}{\left(C_{2E}\textcolor{blue}{\lambda_i} + \textcolor{blue}{h_i} (\sqrt{d_{ij}}-1)\right)^3} $$


2. **Modèle (M2)**

$$Y_{ij} | \psi_i = (\lambda_i,h_i,v_i,pf_i,pv_i)\sim \mathcal{N}(m_{ij},s_{ij}^2)$$

où $$m_{ij} = \frac{1}{E(T_{ij})} \; , \; s_{ij}^2 = \frac{1}{\Delta} \frac{V(T_{ij})}{E(T_{ij})^3}$$


$$E(T_{ij}) = \frac{1}{\textcolor{blue}{pf_i}} \frac{\textcolor{blue}{pv_i}}{1-\textcolor{blue}{pv_i}}\textcolor{blue}{v_i} + C_{2E} \frac{\textcolor{blue}{\lambda_i}}{\textcolor{blue}{pf_i}\sqrt{d_{ij}}-1} + \textcolor{blue}{h_i}$$

$$V(T_{ij}) = \frac{1}{\textcolor{blue}{pf_i}}\left(\frac{\textcolor{blue}{pv_i}}{(1-\textcolor{blue}{pv_i})^2}\right) \textcolor{blue}{v_i^2} + (C_{2V}-C_{2E}^2) \frac{\textcolor{blue}{\lambda_i^2}}{\textcolor{blue}{pf_i}(\sqrt{d_{ij}}-1)^2} + \frac{1-\textcolor{blue}{pf_i}}{\textcolor{blue}{pf_i^2}}\left[\frac{\textcolor{blue}{pv_i}}{1-\textcolor{blue}{pv_i}}\textcolor{blue}{v_i} + C_{2E}\frac{\textcolor{blue}{\lambda_i}}{\sqrt{d_{ij}}-1}\right]^2$$


**versions avec ou sans vigilance (quels paramètres sont concernés?)**

et distributions appropriées sur les paramètres individuels $\psi_i$.

## Estimation

- Implémentation de SAEM sous R pour les deux modèles (M1) et (M2) qui n'entrent pas dans la syntaxe de saemix (du moins, à l'époque...).
- **Observations**
	- Difficulté à estimer certains paramètres, y compris dans le cas de designs riches, en partie à cause des faibles valeurs de $s_{ij}^2$, dues au terme $\Delta$ au dénominateur : motive l'ajout d'un terme résiduel pour modéliser les erreurs de mesure et donne plus de souplesse à l'algorithme dans son exploration. 
	- Pour certaines valeurs de paramètres, on retrouve un comportement classique où les paramètres du modèle sont estimés sans biais ou avec un biais très faible qui diminue à mesure que $N$ (le nombre d'individus observés) augmente, et où la précision des estimations s'améliore à mesure que $N$ augmente. Pour d'autres choix de valeurs de paramètres, certaines estimations sont très biaisées. Une part des explications est peut-être dans la forme des courbes + 2ème partie de l'analyse de sensibilité. **Scenario 1-2-3 à éclaircir**
	- ==(RecapPerdrixMai2019.pdf)== + ==(ResSimusMars.pdf)== + ==(Notes_septembre_2020.pdf)==
- **Quelles valeurs de $\Delta$ a-t-on utilisé dans les simus?**: $\Delta$ est connu, donc non estimé, certaines études sont faites pour $\Delta=10$, (à compléter)...  
- Scripts : 
	- (M1), $V(Y_{ij})= \sigma^2$ ==(SAEM_M1_MD_bruit.R)==
	- (M1), $V(Y_{ij})= s_{ij}^2+\sigma^2$ ==(SAEM_M1_MD_bruit_meca.R)==
	- (M1), $V(Y_{ij})= s_{ij}^2$ ==(SAEM_M1_MD_meca2.R)== 

## Analyse de sensibilité 


Analyse exploratoire pour appréhender si/comment les performances de l'estimation dépendent des valeurs numériques des paramètres du modèle, via différents régimes espérance/variance de la réponse fonctionnelle.

 1. *Quels sont les paramètres qui influencent le plus les sorties du modèle, i.e. la moyenne et la variance de la réponse fonctionnelle, ou le ratio espérance/écart-type de la réponse fonctionnelle ?* ==(Répertoire AS_sorties_modeles)==
		- **Etude graphique** : comparaison graphique de l'allure des courbes de réponses fonctionnelles en fonction des valeurs des paramètres ($H$, $\lambda$,...), des coefficients de variation des effets aléatoires et de $\Delta$ ==(ModelePerdrix_M2avecvigilance.R), (ModelePerdrix_M2avecvigilance_ASeffetsmixtes.R)==  
		- **Etude de sensibilité** avec multisensi, déclinée sous différents modèles 
			1. M1 avec effets aléatoires sur les deux paramètres ==(AS_modeleM1aveceffetsaleatoires.R)==
			2. M2 ==(AS_modeleM2aveceffetsaleatoires.r)==
		
2.  *Quels sont les paramètres qui influencent le plus l'estimation, i.e. le biais et la variance des estimateurs des paramètres ?*
		- Etude menée uniquement dans le modèle M1
		- Pour chaque combinaison de valeurs de paramètres on simule des données et on estime les paramètres (estimateur ponctuel et écart-type - inverse de la matrice d’information de Fisher).
		- Tous les scripts sont dans le sous-répertoire AS_estimation_M1/Pg1 (Pg2 semble être une version préliminaire)
		- Rq: Je ne sais pas quel est le répertoire qui contient les résultats finaux de l'analyse

**Conclusion**: quelques éléments dans ==(PerdrixJanvier2020_reuDynenvie.pdf)==

## Sélection de modèle


1. Dans une version donnée du modèle de réponse fonctionnelle ((M1) ou (M2), version à effets fixes (tous les paramètres sont fixes) ou version à effets aléatoires (tous les paramètres sont aléatoires)), retrouver la bonne forme de la variance parmi:
	- $s_{ij}^2$ (variance mécaniste)
	- $\sigma^2$ (variance résiduelle)
	- $s_{ij}^2+\sigma^2$ (somme des deux)
2. Pour un modèle donné ((M1) ou (M2)), choisir s'il s'agit de sa version à effets fixes ou de sa version à effets mixtes, et avec quelle forme de variance (6 possibilités au total). 


- Quelques éléments sur le design des simulations (valeurs de paramètres, valeurs de N,... ) dans ==(ResM1janv2021.pdf)==