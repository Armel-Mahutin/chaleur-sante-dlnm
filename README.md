# Chaleur et admissions hospitalières : un modèle DLNM reproductible

Ce projet montre, de bout en bout, comment estimer l'effet de la chaleur sur des admissions hospitalières journalières avec un modèle à retards distribués non linéaires (DLNM), puis comment vérifier la solidité du résultat.

Les données sont **simulées**, avec un effet de la chaleur connu à l'avance. Cela permet de vérifier que la méthode retrouve le vrai effet, ce qu'aucune donnée réelle ne permet. La démarche reprend celle de mon stage de recherche à EPIMACT (Inserm U1094 – IRD U270) sur la chaleur et la santé mentale à La Réunion, dont les données ne peuvent pas être diffusées.

Auteur : Armel Mahutin KANHONOU · [Portfolio](https://armel-mahutin.github.io) · ORCID [0009-0005-1764-1327](https://orcid.org/0009-0005-1764-1327)

## Question

Les jours de forte chaleur augmentent-ils les admissions hospitalières, de combien, et sur combien de jours l'effet se prolonge-t-il ?

## Données simulées

- 3 287 jours (2015-2023), environ 27 admissions par jour.
- Température maximale : saisonnalité, persistance d'un jour à l'autre (AR1) et réchauffement de 0,04 °C par an.
- Admissions : saisonnalité propre, effet du jour de la semaine, tendance et surdispersion (loi binomiale négative).
- **Vrai effet** : nul sous 29 °C, puis +3,5 % de risque par °C supplémentaire, réparti sur les jours 0 à 3.

## Méthode

- Régression quasi-Poisson (la surdispersion est estimée, φ ≈ 1,7).
- Base croisée DLNM : spline naturelle à 3 degrés de liberté sur la température, spline naturelle sur les décalages de 0 à 7 jours.
- Ajustement sur la tendance et la saisonnalité (spline du temps, 7 degrés de liberté par an) et sur le jour de la semaine.
- Risque relatif cumulé par rapport à la température médiane (27,1 °C), avec intervalle de confiance à 95 % par la méthode delta.
- La base croisée est programmée avec les seuls packages de base de R (`stats`, `splines`), pour rendre la méthode transparente. Le script `R/99 controle dlnm.R` refait le même modèle avec le package de référence `dlnm`.

## Résultats

| Indicateur | Valeur |
|---|---|
| Risque relatif cumulé à 33,4 °C (P99), estimé | **1,14** (IC 95 % 1,05 – 1,25) |
| Vrai risque relatif simulé | 1,17 |
| Paramètre de surdispersion | 1,69 |

Le modèle retrouve le vrai effet : la valeur simulée est dans l'intervalle de confiance.

![Courbe exposition-réponse](figures/1%20courbe%20exposition%20reponse.png)

![Structure des décalages](figures/2%20structure%20des%20decalages.png)

## Étude de simulation : 200 répétitions

Une seule série ne suffit pas à valider une méthode : le bon résultat peut tenir à la chance. Le script `R/04 etude de simulation.R` simule donc 200 nouvelles séries avec le même mécanisme (graines différentes), refait le même modèle sur chacune et compare le risque relatif cumulé au P99 au vrai effet.

| Indicateur (200 répétitions) | Valeur |
|---|---|
| Vrai RR au P99 (moyenne) | 1,17 |
| RR estimé (moyenne géométrique) | 1,19 |
| Biais relatif du log-RR | +8,9 % |
| RMSE du log-RR | 0,049 |
| Couverture de l'IC à 95 % | **91,5 %** (IC 95 % : 86,7 – 95,0) |
| Puissance (IC excluant RR = 1) | 96,5 % |

**Lecture :** le modèle détecte l'effet de façon fiable (96,5 % des séries), avec une erreur faible. Il le surestime légèrement en moyenne, et son intervalle à 95 % ne contient le vrai effet que dans 91,5 % des cas : il est un peu trop étroit. La cause la plus probable est la forme de l'effet simulé, un seuil net à 29 °C, qu'une spline lisse ne reproduit pas exactement. Avec 4 ou 5 degrés de liberté sur la température, le biais varie entre +6 et +10 % et la couverture entre 92 et 93 % : l'écart persiste. En pratique, cela invite à lire l'intervalle avec prudence et à tester aussi un modèle à seuil.

![Étude de simulation](figures/4%20etude%20de%20simulation.png)

## Sensibilité et diagnostics

Le résultat a été recalculé en faisant varier la durée des décalages (7, 14 et 21 jours) et le contrôle de la saisonnalité (4, 7 et 10 degrés de liberté par an).

- Les 9 estimations vont de 1,10 à 1,24 ; 7 sur 9 excluent l'absence d'effet.
- Allonger les décalages à 14 jours élargit l'intervalle sans gain : la chaleur agit surtout dans les premiers jours.
- Le test de Ljung-Box détecte une faible autocorrélation résiduelle (p ≈ 0,04) ; l'ACF reste très proche de zéro. Une version avec résidu décalé est proposée dans le portfolio.

![Sensibilité et diagnostics](figures/3%20sensibilite%20et%20diagnostics.png)

## Limites

- Données simulées : l'objectif est de montrer la méthode, pas de produire un résultat de santé publique.
- Une seule exposition (température maximale), sans humidité ni pollution de l'air, qui peuvent confondre l'effet de la chaleur.
- Étude écologique : l'association porte sur des jours, pas sur des personnes.

## Reproduire l'analyse

Dans R (4.3 ou plus), depuis le dossier du projet :

```r
source("lancer analyse.R")
```

Packages : `splines` et `stats` (fournis avec R). `dlnm` pour le contrôle croisé facultatif.

## Organisation

```
R/00 fonctions.R           base croisée, risques relatifs cumulés et par décalage
R/01 simuler donnees.R     simulation de la série et du vrai effet
R/02 modele principal.R    modèle principal et figures 1 et 2
R/03 sensibilite.R         9 scénarios de sensibilité, autocorrélation, figure 3
R/04 etude de simulation.R 200 répétitions : biais, RMSE, couverture, figure 4
R/99 controle dlnm.R       même modèle avec le package dlnm
donnees/                   données simulées
resultats/                 tableaux de résultats
figures/                   figures
```

## Référence

Gasparrini A, Armstrong B, Kenward MG. Distributed lag non-linear models. *Statistics in Medicine*. 2010;29(21):2224-2234.
