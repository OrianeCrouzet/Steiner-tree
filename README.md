# Arbre de Steiner dans un graphe, avec et sans restriction

### Choix d'implémentation
Le projet sera codé en Java. 

### Equipe de développement
- AYED Fatima
- CROUZET Oriane

## Graphe géométrique
Un graphe géométrique dans un plan 2D est défini par un ensemble de points dans le plan appelés sommets, et un seuil sur la distance entre les points : il existe une arête entre deux sommets si et seulement si la distance Euclidienne entre les deux sommets est inférieure à ce seuil. Quand une arête existe entre deux sommets, le poids de l’arête est la distance Euclidienne entre les deux sommets dans le plan.

## Problème de l'arbre de Steiner dans un graphe
Etant donnés un graphe G = (VE) et un sous ensemble S inclu dans V de sommets, le problème de l’arbre de Steiner couvrant S consiste à calculer un sous graphe de G qui est un arbre et qui passe par tous les points de S, tel que la longueur totale des arêtes de l’arbre est la plus petite possible.

##  Problème de l’arbre de Steiner avec restriction budgétaire, dans un graphe 
 Soient G = (VE) un graphe, S inclu dans V un sous ensemble de sommets, s de S un point appelé maison-mère, et B un réel appelé budget. Le problème de l’arbre de Steiner couvrant S de budget B passant par s consiste à calculer un sous graphe de G qui est un arbre, qui passe par s, de longueur totale inférieure à B, passant par le plus grand nombre possible de points de S.

## Enoncé du projet
 Il s’agit de proposer une heuristique pour le problème de l’arbre de Steiner dans un graphe géométrique, avec et sans restriction budgétaire. Ce projet s'inscrit dans le cadre de l'UE CPA du master 1 STL de l'université de Sorbonne.
