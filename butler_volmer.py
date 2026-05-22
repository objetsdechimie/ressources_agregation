# -*- coding: utf-8 -*-
"""
Created on Sat Mar 21 12:56:31 2026

@author: maxim
"""

import numpy as np
import matplotlib.pyplot as plt

# Constantes physico-chimiques
R = 8.314  # J/mol/K
F = 96485  # C/mol
T = 298.15 # K

def butler_volmer(eta, j0, alpha, z):
    """Calcule la densité de courant j en fonction de la surtension eta."""
    f = F / (R * T)
    j_anodique = j0 * np.exp(alpha * z * f * eta)
    j_cathodique = -j0 * np.exp(-(1 - alpha) * z * f * eta)
    return j_anodique, j_cathodique, j_anodique + j_cathodique

# Paramètres de simulation
j0 = 1e-3    # A/cm^2 (courant d'échange)
alpha = 0.5  # Facteur de transfert
z = 1        # Nombre d'électrons
surtensions = np.linspace(-0.5, 0.5, 2000) # de -200mV à +200mV

jan, jcat, jtot = butler_volmer(surtensions, j0, alpha, z)

# Tracé
plt.figure(figsize=(8, 5))
plt.plot(surtensions, jtot, 'k', lw=0.05, label='Courant total (j)')
plt.plot(surtensions, jan, '--r', alpha=0.6, label='Composante anodique')
plt.plot(surtensions, jcat, '--b', alpha=0.6, label='Composante cathodique')

plt.axhline(0, color='black', lw=0.5)
plt.axvline(0, color='black', lw=0.5)
plt.xlabel("Surtension $\eta$ (V)")
plt.ylabel("Densité de courant $j$ (A/cm²)")
plt.title("Modèle de Butler-Volmer")
plt.legend()
plt.grid(True, linestyle=':')
plt.show()