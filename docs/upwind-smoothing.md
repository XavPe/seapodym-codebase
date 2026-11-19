# Lissage C1 et C2 du transport près d’une vitesse nulle

Le schéma amont sépare une vitesse signée `v` en deux contributions. Sans
lissage, `v+ = max(v,0)` et `v− = min(v,0)`. Quand `v` change de signe, le
schéma change de branche. Une largeur `epsilon` positive raccorde les deux
branches dans la bande `|v| < epsilon` sans modifier le partage hors bande.
Le calcul direct, les recomputations de l’adjoint et l’adjoint emploient le
même partage. La vitesse est celle utilisée par le solveur de transport.

## Configuration dans un parfile

Placer les éléments directement sous `<par>`, par exemple après `deltaT`.
Les exemples sont des alternatives. Utiliser au plus une déclaration par
élément. Ces réglages numériques sont fixes pendant l’optimisation.

Transport historique, sans lissage :

```xml
<upwind_smoothing_epsilon value="0"/>
```

Lissage quadratique C1 de référence, utilisé dans la campagne RUN007 :

```xml
<upwind_smoothing_mode value="quadratic_C1"/>
<upwind_smoothing_epsilon value="2e-4"/>
```

Lissage quartique C2 à largeur identique au témoin :

```xml
<upwind_smoothing_mode value="quartic_C2"/>
<upwind_smoothing_epsilon value="2e-4"/>
```

C2 avec le même maximum local de diffusion numérique ajoutée que C1 à `2e-4` :

```xml
<upwind_smoothing_mode value="quartic_C2"/>
<upwind_smoothing_epsilon value="2.6666666666666667e-4"/>
```

En l’absence de `upwind_smoothing_mode`, la station utilise `quadratic_C1`.
La largeur est lue uniquement dans le parfile XML. En l’absence de
`upwind_smoothing_epsilon`, elle vaut zéro et le lissage de vitesse
est désactivé. Les variables d’environnement sont ignorées.
Déclarer la largeur et le mode dans le parfile pour reproduire un calcul.
Une largeur négative, non finie ou illisible, ou un mode inconnu, arrête
l’initialisation. Le journal indique le mode, la largeur et son origine.
Les réglages effectifs figurent aussi dans le JSON du mode local `--fd-point`.

## Définition et effet sur la diffusion

Soit `r = v/epsilon` dans la bande. Les deux modes utilisent
`v+ = (v + phi)/2` et `v− = (v - phi)/2`, avec `v+ + v− = v`.

| Mode | `phi` dans la bande | Raccord à `|v| = epsilon` | `phi` à `v = 0` |
|---|---|---|---|
| `quadratic_C1` | `epsilon * (1 + r²) / 2` | valeur et dérivée première continues | `epsilon / 2` |
| `quartic_C2` | `epsilon * (3 + 6*r² - r⁴) / 8` | valeur, dérivées première et seconde continues | `3*epsilon / 8` |

Dans les deux cas, `phi = |v|` hors bande et les branches arithmétiques
historiques sont conservées. À largeur égale, le maximum local de diffusion
numérique supplémentaire de C2 vaut les trois quarts de celui de C1.
Le facteur `4/3` sur la largeur de C2 égalise ces deux maxima. Il change
aussi la largeur de la bande et ne reproduit pas le profil de diffusion C1.

Les valeurs ci-dessus sont des points de comparaison testés, à réévaluer
avec la grille et le pas de temps si la configuration change. Dans la
campagne ERA5 depuis RUN007, aucun gain de C2 sur C1 n'a été établi.
La régularité du raccord local ne garantit pas celle de tous les termes
du modèle ou de l’objectif. Pour l’option qui adoucit le gradient d’habitat
dans la diffusion, voir [Régularisation du gradient d’habitat](habitat-gradient-smoothing.md).
