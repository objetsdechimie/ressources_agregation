# -*- coding: utf-8 -*-
"""
Created on Sat Mar 21 13:04:06 2026

@author: maxim
"""

import numpy as np
import matplotlib.pyplot as plt

def solve_huckel(matrix, name):
    """Calcule et trie les niveaux d'énergie d'une matrice de Hückel."""
    energies = np.linalg.eigvalsh(matrix)
    # On trie du plus stable (énergie la plus basse) au moins stable
    # Rappel : E = alpha + x*beta. Comme beta < 0, le x le plus grand est le plus stable.
    return sorted(energies, reverse=True)

# 1. Éthène (2 carbones, 1 liaison double)
H_ethene = np.array([[0, 1],
                     [1, 0]])

# 2. Benzène (Cycle à 6 carbones)
H_benzene = np.zeros((6, 6))
for i in range(6):
    H_benzene[i, (i+1)%6] = 1
    H_benzene[(i+1)%6, i] = 1

# 3. Butadiène (Chaîne linéaire à 4 carbones)
H_butadiene = np.zeros((4, 4))
for i in range(3):
    H_butadiene[i, i+1] = 1
    H_butadiene[i+1, i] = 1

# Calculs
molecules = {
    "Éthène": solve_huckel(H_ethene, "Éthène"),
    "Benzène": solve_huckel(H_benzene, "Benzène"),
    "Butadiène": solve_huckel(H_butadiene, "Butadiène")
}

# Visualisation
fig, axs = plt.subplots(1, 3, figsize=(12, 6), sharey=True)
fig.suptitle("Niveaux d'énergie de Hückel ($E = \\alpha + x\\beta$)", fontsize=14)

for ax, (name, levels) in zip(axs, molecules.items()):
    # On trace chaque niveau comme une ligne horizontale
    for e in levels:
        ax.hlines(e, 0.2, 0.8, colors='blue', lw=2)
        ax.text(0.85, e, f"{e:.2f}", va='center', fontsize=9)
    
    ax.set_title(name)
    ax.set_xticks([])
    ax.set_ylabel("Coefficient $x$ de $\\beta$")
    ax.grid(axis='y', linestyle='--', alpha=0.7)
    # Limites pour bien voir les niveaux
    ax.set_ylim(-2.5, 2.5)

plt.tight_layout()
plt.show()