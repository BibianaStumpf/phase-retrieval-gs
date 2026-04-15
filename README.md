# Récupération de Phase — Algorithme de Gerchberg–Saxton

> Simulation numérique de la récupération de phase en optique par l'algorithme itératif de Gerchberg–Saxton.

---

## Contexte physique

En optique, un champ électromagnétique est une grandeur complexe :

$$E(r) = A(r)\,e^{i\phi(r)}$$

où $A(r)$ est l'amplitude (liée à l'intensité) et $\phi(r)$ la phase (liée au déphasage du front d'onde). Les capteurs optiques (CCD, CMOS…) ne mesurent que l'**intensité** :

$$I(r) = |E(r)|^2$$

Cette mesure entraîne une **perte totale d'information sur la phase**, rendant impossible la reconstruction directe du champ incident. C'est le **problème de récupération de phase**, fondamental en astronomie, cristallographie, microscopie et profilométrie.

---

## Algorithme de Gerchberg–Saxton

L'algorithme G–S est une méthode itérative qui alterne entre le plan objet et le plan de Fourier, en imposant à chaque étape les contraintes physiques connues :

```
x₀ = |f| · exp(iφ₀)          ← initialisation (phase aléatoire ou nulle)
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

## Contenu du dépôt

```
phase-retrieval-GS/
├── gerchberg_saxton.ipynb   # Notebook principal (implémentation + expériences)
└── README.md
```

---

## Expériences réalisées

### 1. Phases discontinues

L'objet est un disque de 120×120 pixels avec deux zones de phase distinctes (anneau et centre). On augmente progressivement le saut de phase pour caractériser les limites de l'algorithme.

| Cas | Phase anneau | Phase centre | Saut total | Résultat |
|-----|-------------|-------------|------------|---------|
| 1 | 0 | 0 | 0 | Reconstruction parfaite |
| 2 | π/2 | −π/2 | π | Structure correctement retrouvée |
| 3 | π | −π/2 | 3π/2 | Début de dégradation |
| 4 | 7π/4 | 7π/4 | — | Ambiguïté mod 2π illustrée |

### 2. Rampe de phase (onde plane inclinée)

Une onde plane inclinée de vecteur $(k_x, k_y)$ introduit une rampe de phase spatiale. On vérifie que cela décale le pic d'intensité dans le plan de Fourier, et on teste la capacité de reconstruction de l'algorithme.

| Cas | $k_x$ | $k_y$ | Résultat |
|-----|--------|--------|---------|
| Référence | 0 | 0 | Phase uniforme |
| Rampe horizontale | 0.1 | 0 | Bandes verticales reconstruites |
| Rampe diagonale | 0.3 | 0.3 | Décalage diagonal confirmé |
| Rampe + discontinuité | 1.0 | 0.5 | Aliasing + stagnation |

---

## Limites identifiées

| Phénomène | Origine |
|-----------|---------|
| Bruit dans les zones d'amplitude nulle | Limite numérique (arg de nombres ≈ 0) |
| Ambiguïté modulo 2π | Invariance de l'intensité |
| Perte de phase absolue | Invariance par translation de phase globale |
| Stagnation (minima locaux) | Nature non-convexe du problème |
| Repliement de phase (aliasing) | Variation > π entre pixels adjacents |
| Vortex de phase | Singularité physique (phase indéfinie si amplitude = 0) |

---

## Technologies

- **Python** · NumPy · Matplotlib
- **Transformée de Fourier rapide (FFT)** : `numpy.fft`
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

*Projet réalisé dans le cadre du cours d'optique — École Centrale Méditerranée, 2025*
