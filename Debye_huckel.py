# -*- coding: utf-8 -*-
"""
Created on Sat Mar 21 12:58:48 2026

@author: maxim
"""

import numpy as np
import matplotlib.pyplot as plt

def debye_huckel(I, z, a=3e-10):
    """
    Calcule gamma pour un ion de charge z.
    a : paramètre de taille de l'ion (m) - par défaut 3 Angstroms
    """
    A = 0.509  # L^1/2.mol^-1/2 à 25°C
    B = 3.29e9 # m^-1.L^1/2.mol^-1/2 à 25°C
    
    log_gamma = - (A * z**2 * np.sqrt(I)) / (1 + B * a * np.sqrt(I))
    return 10**log_gamma

# Simulation pour différents ions
force_ionique = np.linspace(0.0001, 0.5, 200)

gamma_z1 = debye_huckel(force_ionique, z=1)
gamma_z2 = debye_huckel(force_ionique, z=2)
gamma_z3 = debye_huckel(force_ionique, z=3)

# Tracé
plt.figure(figsize=(8, 5))
plt.plot(force_ionique, gamma_z1, label='Ion charge |z|=1 (ex: Na+)')
plt.plot(force_ionique, gamma_z2, label='Ion charge |z|=2 (ex: Mg2+)')
plt.plot(force_ionique, gamma_z3, label='Ion charge |z|=3 (ex: Al3+)')

plt.xlabel("Force ionique $I$ (mol/L)")
plt.ylabel("Coefficient d'activité $\gamma_i$")
plt.title("Loi de Debye-Hückel étendue")
plt.legend()
plt.grid(True)
plt.show()