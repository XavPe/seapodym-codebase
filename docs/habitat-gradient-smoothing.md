# Régularisation du gradient d’habitat dans la diffusion

Lorsque la qualité de l’habitat change d’une maille à la suivante, le modèle
réduit la diffusion de façon à tenir compte de l’amplitude du changement.
La valeur absolue employée dans ce calcul présente une pointe lorsque le
gradient passe par zéro. Cette option remplace la pointe par un raccord
lisse dans le calcul de diffusion des adultes.

Pour chaque direction, soit `g` le gradient d’habitat et `d` le pas de grille
effectif dans cette direction. Le produit `q = d*g` est exprimé dans les
unités de l’indice d’habitat. Le modèle remplace `|q|` par
`sqrt(q² + eta²) - eta` lorsque `eta > 0`, avec une forme de calcul équivalente
pour éviter de soustraire deux nombres presque égaux près de zéro. La valeur
reste nulle pour `q = 0` et varie sans angle à son voisinage. La valeur
adoucie est inférieure ou égale à `|q|`, ce qui atténue localement la
réduction de diffusion. L’adjoint applique la dérivée correspondante.
Le gradient signé qui dirige la nage garde sa formule. Le champ d’habitat
et les paramètres de vraisemblance gardent leurs définitions.

## Configuration dans un parfile

Ajouter directement sous `<par>`, par exemple après `deltaT` :

```xml
<habitat_gradient_smoothing_eta value="1e-7"/>
```

`1e-7` est la valeur positive employée dans les optimisations de comparaison
depuis RUN007. `1e-8` et `1e-6` ont servi de points bas et haut dans les
comparaisons à paramètres fixes. `eta` est une largeur en unités de `q`,
indépendante de la largeur de vitesse `upwind_smoothing_epsilon`. Utiliser
une valeur finie et non négative. L’absence de l’élément ou `value="0"`
conserve la modulation de diffusion historique, avec sa branche adjointe.

Exemple combiné, C1 de référence et habitat lissé :

```xml
<upwind_smoothing_mode value="quadratic_C1"/>
<upwind_smoothing_epsilon value="2e-4"/>
<habitat_gradient_smoothing_eta value="1e-7"/>
```

Pour isoler l’habitat, conserver le même réglage de vitesse dans les deux
parfiles et comparer `eta=0` à `eta=1e-7` au même vecteur de paramètres.
Le réglage est fixe et la valeur effective figure au démarrage ainsi que
dans le JSON du mode local `--fd-point`. Aucun gain sur le témoin n’a été
établi dans la campagne étudiée. Ces valeurs restent des exemples de
sensibilité. Voir aussi [Lissage de vitesse C1 et C2](upwind-smoothing.md).
