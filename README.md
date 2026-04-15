# Récupération de Phase — Algorithme de Gerchberg–Saxton

> Simulation numérique de la récupération de phase en optique par l'algorithme itératif de Gerchberg–Saxton.

---

## Contexte physique

En optique, un champ électromagnétique est une grandeur complexe :

$$E(r) = A(r)e^{i\phi(r)}$$

Les capteurs (CCD, CMOS…) ne mesurent que l'intensité :

$$I(r) = |E(r)|^2$$

Cette mesure entraîne une **perte totale d'information sur la phase**. C'est le **problème de récupération de phase**, fondamental en astronomie, cristallographie, microscopie et profilométrie.

---

## Algorithme de Gerchberg–Saxton

```
x₀ = |f| · exp(iφ₀)          ← initialisation (phase nulle ou aléatoire)
  ↓  FFT
Xₖ = F(xₖ)
  ↓  contrainte Fourier
X'ₖ = |F| · exp(iψₖ)          ← module mesuré, phase calculée conservée
  ↓  IFFT
x'ₖ = F⁻¹(X'ₖ)
  ↓  contrainte objet
xₖ₊₁ = |f| · exp(iφₖ)        ← amplitude connue, phase estimée conservée
```

---

## Structure du dépôt

```
phase-retrieval-GS/
├── README.md
├── notebooks/
│   └── gerchberg_saxton.ipynb   ← notebook principal
├── src/
│   └── gs_algorithm.py          ← fonctions réutilisables
└── results/
    └── figures/                 ← figures générées automatiquement
```

---

## Contenu du notebook

| Section | Description |
|---------|-------------|
| 0. Configuration | Imports, paramètres globaux, création du dossier `results/figures/` |
| 1. Objet de référence | Disque binaire avec deux zones de phase |
| 2. Phase initiale | Comparaison convergence : phase nulle vs. aléatoire |
| 3. Phases discontinues | Sauts croissants, ambiguïté mod 2π |
| 4. Phase gaussienne | Profil radial continu — cas plus réaliste |
| 5. Rampe de phase | Onde plane inclinée, décalage dans le plan de Fourier |
| 6. Tableau de synthèse | Pearson r et RMSE pour tous les cas |
| 7. Analyse des limites | Illustrations des 6 limites fondamentales |
| 8. Conclusion | Comparaison avec HIO, DPR/MPR |

---

## Métriques de qualité

Deux métriques sont calculées automatiquement pour chaque reconstruction, **uniquement dans le masque** (zone d'amplitude non nulle) :

- **Corrélation de Pearson** $r$ : mesure la ressemblance structurelle entre phase réelle et estimée (1 = parfait)
- **RMSE** (rad) : erreur quadratique moyenne de phase après alignement global

```python
from src.gs_algorithm import pearson_phase, rmse_phase
r   = pearson_phase(phase_reelle, phase_estimee, masque=amplitude > 0)
rms = rmse_phase(phase_reelle, phase_estimee, masque=amplitude > 0)
```

---

## Limites identifiées

| Phénomène | Origine | Contournable ? |
|-----------|---------|---------------|
| Bruit dans les zones nulles | Limite numérique (arg ≈ 0) | Oui — masque |
| Ambiguïté mod 2π | Invariance de l'intensité | Non (physique) |
| Perte de phase absolue | Invariance par translation globale | Non (physique) |
| Stagnation (minima locaux) | Problème non-convexe | Partiellement — HIO |
| Repliement de phase (aliasing) | Variation > π entre pixels | Partiellement — sur-échantillonnage |
| Vortex de phase | Singularité (amplitude = 0) | Non (physique) |

---

## Technologies

- **Python** · NumPy · Matplotlib
- **FFT** : `numpy.fft`
- Environnement : Google Colab / Jupyter Notebook

---

## Références

- Fienup, J. R. (1982). *Phase retrieval algorithms: a comparison.* Applied Optics, 21(15), 2758–2769.
- Zhao, T. & Chi, Y. (2020). *Modified Gerchberg–Saxton (G-S) Algorithm and Its Application.* Entropy, 22(12), 1354.
- Thimons, T. & Wittle, L. (2018). *Investigating the Gerchberg-Saxton Phase Recovery Algorithm.*
- Carpenter, J. *Gerchberg Saxton algorithm (Tutorial).* The University of Queensland.

---

## Auteur

**Bibiana Terres Stumpf**  
Double diplôme — Ingénierie Généraliste (École Centrale Méditerranée) & Ingénierie Physique (UFRGS)

*Projet réalisé dans le cadre du cours d'optique avec ALSENE Adrien, LEMONSU Hugo, DOMMELIER Maylis et STEFANOVIC Katarina -— École Centrale Méditerranée, 2025*
