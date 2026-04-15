"""
gs_algorithm.py
===============
Fonctions principales pour la récupération de phase par l'algorithme de Gerchberg–Saxton.

Utilisation depuis le notebook :
    import sys
    sys.path.append('../src')
    from gs_algorithm import gerchberg_saxton, construire_objet, pearson_phase
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import os


def make_grid(N):
    """Retourne les grilles X, Y, R pour une image N×N centrée en (0,0)."""
    x = np.arange(N) - N / 2
    X, Y = np.meshgrid(x, x)
    R = np.sqrt(X**2 + Y**2)
    return X, Y, R


def construire_objet(N, rayon_ext, rayon_int, phase_anneau, phase_centre):
    """
    Construit un champ complexe avec deux zones de phase (anneau + centre).

    Paramètres
    ----------
    N : int — taille de la grille (N×N)
    rayon_ext : float — rayon extérieur du disque
    rayon_int : float — rayon intérieur (centre)
    phase_anneau : float — phase de l'anneau (rad)
    phase_centre : float — phase du centre (rad)

    Retourne
    --------
    amplitude : ndarray (N×N) — amplitude binaire
    champ     : ndarray (N×N, complexe) — champ avec phase imposée
    """
    _, _, R = make_grid(N)

    amplitude = np.zeros((N, N))
    amplitude[R < rayon_ext] = 1.0

    phase = np.zeros((N, N))
    phase[R < rayon_ext] = phase_anneau
    phase[R < rayon_int] = phase_centre

    return amplitude, amplitude * np.exp(1j * phase)


def construire_objet_gaussien(N, rayon_ext, A=1.0, r0=0.0, sigma=10.0):
    """
    Construit un champ dont la phase suit un profil gaussien radial.

    Phase : phi(r) = A * exp(-((r - r0)^2) / sigma^2)

    Paramètres
    ----------
    N        : int   — taille de la grille
    rayon_ext: float — rayon extérieur du disque
    A        : float — amplitude de la gaussienne (rad)
    r0       : float — centre radial de la gaussienne
    sigma    : float — largeur de la gaussienne

    Retourne
    --------
    amplitude, champ
    """
    _, _, R = make_grid(N)

    amplitude = np.zeros((N, N))
    amplitude[R < rayon_ext] = 1.0

    phase = A * np.exp(-((R - r0)**2) / sigma**2)
    phase[R >= rayon_ext] = 0.0

    return amplitude, amplitude * np.exp(1j * phase)


def construire_objet_rampe(N, rayon_ext, rayon_int, phase_anneau, phase_centre, kx=0.0, ky=0.0):
    """
    Construit un champ avec une rampe de phase (onde plane inclinée).

    E(x,y) = A(x,y) * exp(i*(phase_zones + kx*x + ky*y))

    Paramètres
    ----------
    kx, ky : float — composantes du vecteur d'onde (contrôlent l'inclinaison)

    Retourne
    --------
    amplitude, champ
    """
    X, Y, R = make_grid(N)

    amplitude = np.zeros((N, N))
    amplitude[R < rayon_ext] = 1.0

    phase = np.zeros((N, N))
    phase[R < rayon_ext] = phase_anneau
    phase[R < rayon_int] = phase_centre

    rampe = np.exp(1j * (kx * X + ky * Y))
    champ = amplitude * np.exp(1j * phase) * rampe
    champ[R >= rayon_ext] = 0.0

    return amplitude, champ


# ─────────────────────────────────────────────
# Algorithme de Gerchberg–Saxton
# ─────────────────────────────────────────────

def gerchberg_saxton(amplitude_objet, TF_mesure, n_iter=500, phase_init='zeros'):
    """
    Algorithme de Gerchberg–Saxton pour la récupération de phase.

    Itère entre le plan objet et le plan de Fourier en imposant à chaque
    étape les contraintes physiques connues.

    Paramètres
    ----------
    amplitude_objet : ndarray (N×N, réel)
        Amplitude connue dans le plan objet.
    TF_mesure : ndarray (N×N, réel)
        Module de la TF mesuré (= sqrt(intensité) expérimentale).
    n_iter : int
        Nombre d'itérations.
    phase_init : str
        'zeros'  → phase initiale nulle
        'random' → phase initiale aléatoire uniforme dans [0, 2π)

    Retourne
    --------
    champ : ndarray (N×N, complexe)
        Champ estimé dans le plan objet après convergence.
    erreurs_fourier : list[float]
        Erreur quadratique normalisée dans le plan de Fourier.
    erreurs_objet : list[float]
        Erreur quadratique dans le plan objet.
    """
    N = amplitude_objet.shape[0]

    if phase_init == 'random':
        phi0 = 2 * np.pi * np.random.rand(N, N)
    else:
        phi0 = np.zeros((N, N))

    champ = amplitude_objet * np.exp(1j * phi0)
    erreurs_fourier = []
    erreurs_objet   = []

    for _ in range(n_iter):
        # 1. TF directe
        G = np.fft.fftshift(np.fft.fft2(champ))

        # Erreur Fourier (normalisée)
        erreurs_fourier.append(np.sqrt(np.sum((np.abs(G) - TF_mesure)**2) / N**2))

        # 2. Contrainte Fourier : imposer le module mesuré
        G2 = TF_mesure * np.exp(1j * np.angle(G))

        # 3. TF inverse
        g = np.fft.ifft2(np.fft.ifftshift(G2))

        # Erreur objet
        erreurs_objet.append(np.sqrt(np.sum((np.abs(g) - amplitude_objet)**2)))

        # 4. Contrainte objet : imposer l'amplitude connue
        champ = amplitude_objet * np.exp(1j * np.angle(g))

    return champ, erreurs_fourier, erreurs_objet


# ─────────────────────────────────────────────
# Métriques de qualité
# ─────────────────────────────────────────────

def pearson_phase(phase_reelle, phase_estimee, masque=None):
    """
    Corrélation de Pearson entre la phase réelle et estimée.

    On ne calcule la corrélation que dans les zones d'amplitude non nulle
    (masque). La phase est définie modulo 2π donc on compare
    exp(i*phi) pour être insensible aux constantes globales.

    Paramètres
    ----------
    phase_reelle   : ndarray (N×N) — phase de référence (rad)
    phase_estimee  : ndarray (N×N) — phase reconstruite (rad)
    masque         : ndarray (N×N, bool) ou None — zone d'intérêt

    Retourne
    --------
    r : float — coefficient de Pearson dans [-1, 1]
    """
    if masque is not None:
        p_r = phase_reelle[masque].ravel()
        p_e = phase_estimee[masque].ravel()
    else:
        p_r = phase_reelle.ravel()
        p_e = phase_estimee.ravel()

    # Alignement global : minimiser l'erreur par rapport à un décalage constant
    delta = np.angle(np.sum(np.exp(1j * (p_r - p_e))))
    p_e_alignee = p_e + delta

    r = np.corrcoef(p_r, p_e_alignee)[0, 1]
    return float(r)


def rmse_phase(phase_reelle, phase_estimee, masque=None):
    """
    RMSE de l'erreur de phase (en radians), après alignement global.

    Paramètres
    ----------
    phase_reelle, phase_estimee : ndarray (N×N)
    masque : ndarray (N×N, bool) ou None

    Retourne
    --------
    rmse : float (rad)
    """
    if masque is not None:
        p_r = phase_reelle[masque].ravel()
        p_e = phase_estimee[masque].ravel()
    else:
        p_r = phase_reelle.ravel()
        p_e = phase_estimee.ravel()

    delta = np.angle(np.sum(np.exp(1j * (p_r - p_e))))
    err   = np.angle(np.exp(1j * (p_r - (p_e + delta))))
    return float(np.sqrt(np.mean(err**2)))


# ─────────────────────────────────────────────
# Visualisation
# ─────────────────────────────────────────────

def afficher_resultats(phase_reelle_2D, champ_estime, erreurs_fourier, erreurs_objet,
                       titre='', masque=None, save_path=None):
    """
    Affiche côte à côte : phase réelle, phase estimée, convergence, coupe centrale.
    Calcule et affiche les métriques de Pearson et RMSE.

    Paramètres
    ----------
    save_path : str ou None — si fourni, sauvegarde la figure à ce chemin
    """
    phase_estimee_2D = np.angle(champ_estime)
    N = phase_reelle_2D.shape[0]

    r   = pearson_phase(phase_reelle_2D, phase_estimee_2D, masque)
    rms = rmse_phase(phase_reelle_2D, phase_estimee_2D, masque)

    fig = plt.figure(figsize=(18, 4))
    full_title = f"{titre}   |   Pearson r = {r:.4f}   |   RMSE = {rms:.4f} rad"
    fig.suptitle(full_title, fontsize=12, y=1.01)
    gs_layout = gridspec.GridSpec(1, 4, figure=fig)

    # Phase réelle
    ax0 = fig.add_subplot(gs_layout[0])
    im0 = ax0.imshow(phase_reelle_2D, cmap='hsv', vmin=-np.pi, vmax=np.pi)
    ax0.set_title('Phase réelle')
    ax0.axis('off')
    _colorbar_phase(fig, im0, ax0)

    # Phase estimée
    ax1 = fig.add_subplot(gs_layout[1])
    im1 = ax1.imshow(phase_estimee_2D, cmap='hsv', vmin=-np.pi, vmax=np.pi)
    ax1.set_title('Phase estimée (G–S)')
    ax1.axis('off')
    _colorbar_phase(fig, im1, ax1)

    # Convergence
    ax2 = fig.add_subplot(gs_layout[2])
    ax2.semilogy(erreurs_fourier, label='Plan Fourier', linewidth=1.5)
    ax2.semilogy(erreurs_objet,   label='Plan objet',   linewidth=1.5)
    ax2.set_xlabel('Itération')
    ax2.set_ylabel('Erreur (log)')
    ax2.set_title('Convergence')
    ax2.legend(fontsize=8)
    ax2.grid(True, which='both', alpha=0.4)

    # Coupe centrale
    ax3 = fig.add_subplot(gs_layout[3])
    cx  = np.linspace(-0.5, 0.5, N)
    row = N // 2
    ax3.plot(cx, phase_reelle_2D[row, :],  'b.', markersize=3, label='Réelle')
    ax3.plot(cx, phase_estimee_2D[row, :], 'r.', markersize=3, label='Estimée')
    ax3.set_ylim(-np.pi - 0.5, np.pi + 0.5)
    ax3.set_xlabel('Position normalisée')
    ax3.set_ylabel('Phase (rad)')
    ax3.set_title('Coupe centrale')
    ax3.legend(fontsize=8)
    ax3.grid(True, alpha=0.4)

    plt.tight_layout()

    if save_path:
        os.makedirs(os.path.dirname(save_path), exist_ok=True)
        plt.savefig(save_path, dpi=150, bbox_inches='tight')

    plt.show()
    return r, rms


def _colorbar_phase(fig, im, ax):
    cbar = fig.colorbar(im, ax=ax, shrink=0.8)
    cbar.set_ticks([-np.pi, 0, np.pi])
    cbar.set_ticklabels([r'$-\pi$', '0', r'$\pi$'])
    cbar.set_label('Phase (rad)', fontsize=8)
