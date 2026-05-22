# -*- coding: utf-8 -*-
"""
Created on Thu Mar 12 11:15:18 2026

@author: maxim
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider

alpha = -8.1  # Energie arbitraire de l'orbitale consideree
beta = -1     # Valeur de l'integrale de recouvrement entre deux orbitales voisines (arbitraire)
n_max = 100000
a = 1e-10
kb = 1.35e-23 #J/K
T = 298 #K

def E(N, n):
    # Renvoie l'energie du n-ieme niveau d'energie d'une chaine de N atomes
    return alpha + 2 * beta * np.cos(n * np.pi / (N + 1))

def compute_energies(N):
    return np.array([E(N, n) for n in range(1, N + 1)])

Emin = alpha + 2 * beta   # = alpha - 2|beta| (bas de bande)
Emax = alpha - 2 * beta   # = alpha + 2|beta| (haut de bande)

def dos_analytique(e_range, N):
    """
    Densité d'états analytique exacte pour une chaîne de Hückel 1D (limite N→∞) :
        g(ε) = N/π · 1/sqrt(4β² - (ε-α)²)
    C'est une fonction en U avec divergences de Van Hove aux deux bords.
    On clip l'intérieur du sqrt pour éviter la division par zéro.
    """
    inside = 4 * beta**2 - (e_range - alpha)**2
    inside = np.clip(inside, 1e-6, None) #evite une idvision par zero 
    return N / np.pi / np.sqrt(inside)



# --- 3 Subplots ---
fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(18, 6)) #définit nos 3 figures
plt.subplots_adjust(bottom=0.25, wspace=0.4)

# -------------------------------------------------------
# AX1 : Niveaux d'énergie (code original)
# -------------------------------------------------------
A = [0] * n_max #créé liste de 1000 zeros 
N_init = 1

for i in range(0, n_max):
    if i < N_init:
        Enew = E(N_init, i + 1)
        if i == 0:
            A[i], = ax1.plot([-i, i], [Enew, Enew], color='black', label='occupé') #trait horizontal centré en 0 dont largeur croit avec indice i
        else:
            A[i], = ax1.plot([-i, i], [Enew, Enew], color='white')
    else:
        if i == N_init:
            A[i], = ax1.plot([-i, i], [Enew, Enew ], color='red', label='libre') 
        else:
            A[i], = ax1.plot([-i, i], [Enew , Enew ], color='white')

ax1.set_ylim(-10.4, -5.8)
ax1.set_xlabel('k', fontsize=15)

ax1.set_ylabel('Energie (eV)', fontsize=15)

ax1.set_title('Niveaux électroniques (Hückel)', fontsize=13)
ax1.tick_params(axis='x', colors='white')
ax1.tick_params(axis='y', labelsize=13)
ax1.legend(loc='upper right')

# -------------------------------------------------------
# AX2 : Histogramme + courbe DOS analytique
# -------------------------------------------------------
def draw_histogram(N):
    ax2.cla() #clear the axes
    energies = compute_energies(N) #calcules les niveaux expérimentales
    n_bins = min(N, 100) #entier
    ax2.hist(energies, bins=n_bins, color='steelblue', edgecolor='black',
             alpha=0.6, orientation='vertical', density=False)

    if N >= 10:
        margin = 0.005 * abs(2 * beta) #pour pas coler au bord
        e_range = np.linspace(Emin + margin, Emax - margin, 2000)
        g = dos_analytique(e_range, N)
        bin_width = (Emax - Emin) / n_bins #largeur d'une bande
        # La courbe analytique * bin_width donne le nombre de niveaux par bin
        ax2.plot(e_range, g * bin_width, color='red', lw=2.5, label='DOS analytique')
        ax2.legend(fontsize=11)

    ax2.set_xlabel('Energie (eV)', fontsize=15)
    ax2.set_ylabel('Nombre de niveaux', fontsize=15)
    ax2.set_title("Densité des niveaux d'énergie", fontsize=13)
    ax2.set_xlim(Emin - 0.1, Emax + 0.1)
    ax2.tick_params(labelsize=13)

draw_histogram(N_init)

# -------------------------------------------------------
# AX3 : g(ε) depuis histogramme lissé, ε en X, g en Y
# -------------------------------------------------------
def draw_dos_plot(N, ef_idx):
    ax3.cla()

    energies = compute_energies(N)
    n_bins = min(N, 100) #decoupe le plot en intervalle egaux : cb d eniveaux d'energies dans un intervalle

    counts, bin_edges = np.histogram(energies, bins=n_bins)
    bin_centers = 0.5 * (bin_edges[:-1] + bin_edges[1:])

    n_occupied = min(ef_idx, N)
    E_fermi = E(N, n_occupied) if n_occupied > 0 else Emin

    if N >= 10:
        from scipy.ndimage import uniform_filter1d
        smooth_window = max(3, n_bins // 15) #taille fenetre de lissage 
        counts_smooth = uniform_filter1d(counts.astype(float), size=smooth_window) #uniform permet de filtrer
        #pour chaque bin, on remplace sa valeur par la moyenne des `smooth_window` bins voisins.
        #Cela évite que la courbe dans AX3 soit trop irrégulière/crénelée quand N est petit.

        # Interpolation dense
        e_dense = np.linspace(bin_centers[0], bin_centers[-1], 1000)
        g_dense = np.interp(e_dense, bin_centers, counts_smooth)

        # g en X, ε en Y
        ax3.plot(g_dense, e_dense, color='black', lw=2.5)

        # Hachures à gauche de εF
        mask = e_dense <= E_fermi
        ax3.fill_betweenx(e_dense, g_dense, where=mask,
                          facecolor='none', edgecolor='red',
                          linewidth=0.8, hatch='////', label='états occupés')

        # Ligne horizontale εF
        ax3.axhline(y=E_fermi, color='red', lw=1.5, linestyle='--') #trace ligne horizontale 
        ax3.text(max(g_dense) * 1.02, E_fermi, r'$\varepsilon_F$',
                 color='red', fontsize=14, ha='left', va='center', clip_on=False)

    elif N >= 2:
        ax3.plot(bin_centers, counts, color='black', lw=2.5)
    else:
        ax3.axvline(x=energies[0], color='steelblue', lw=1.5)

    ax3.set_xlabel(r'$g(\varepsilon)$', fontsize=15)
    ax3.set_ylabel(r'$\varepsilon$ (eV)', fontsize=15)
    ax3.set_title(r"Densité d'états $g(\varepsilon)$", fontsize=13)
    ax3.set_ylim(Emin - 0.1, Emax + 0.1)
    ax3.set_xlim(left=0)
    ax3.tick_params(labelsize=13)
    ax3.legend(fontsize=11)

draw_dos_plot(N_init, ef_idx=N_init // 2 + 1)

# -------------------------------------------------------
# Slider
# -------------------------------------------------------
axN = plt.axes([0.25, 0.1, 0.65, 0.03])
N_slider = Slider(
    ax=axN,
    label="Nombres d'atomes",
    valmin=1,
    valmax=n_max,
    valstep=1,
    valinit=N_init)

def fermi(En, ef, T, kb):
    return 1 / (1 + np.exp((En - ef) / (kb * T)))

def update(val):
    N = N_slider.val
    ef = int(val) // 2 + 1

    # Mise à jour ax1
    for i in range(0, n_max):
        if ef < i:
            if i < N:
                Enew = E(N, i + 1)
                A[i].set_ydata([Enew, Enew])
                A[i].set_color('red')
            if i >= N:
                A[i].set_ydata([0, 0])
                A[i].set_color('white')
        else:
            if i < N:
                Enew = E(N, i + 1)
                A[i].set_ydata([Enew, Enew])
                A[i].set_color('black')
            if i >= N:
                A[i].set_ydata([0, 0])
                A[i].set_color('white')

    draw_histogram(N)
    draw_dos_plot(N, ef_idx=ef)
    fig.canvas.draw_idle()

N_slider.on_changed(update)

fig.suptitle('Energie des niveaux électroniques (méthode de Hückel)', fontsize=16, y=0.98)
plt.show()