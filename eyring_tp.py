# -*- coding: utf-8 -*-
"""
Created on Sat Mar 21 12:58:22 2026

@author: maxim
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import linregress

# Données expérimentales (Exemple : k en s^-1, T en K)
T_exp = np.array([283, 293, 303, 313, 323])
k_exp = np.array([1.2e-4, 4.5e-4, 1.5e-3, 4.8e-3, 1.4e-2])

# Constantes
R = 8.314
kB = 1.38e-23
h = 6.626e-34

# Variables pour la régression
x = 1 / T_exp
y = np.log(k_exp / T_exp)

# Régression linéaire
slope, intercept, r_value, p_value, std_err = linregress(x, y)

# Calcul des grandeurs d'activation
dH_activ = -slope * R
dS_activ = (intercept - np.log(kB / h)) * R

print(f"Enthalpie d'activation : {dH_activ/1000:.2f} kJ/mol")
print(f"Entropie d'activation  : {dS_activ:.2f} J/mol/K")
print(f"Coefficient de corrélation R² : {r_value**2:.4f}")

# Visualisation
plt.scatter(x, y, color='red', label='Données exp.')
plt.plot(x, slope * x + intercept, label='Régression d\'Eyring')
plt.xlabel("1/T ($K^{-1}$)")
plt.ylabel("ln(k/T)")
plt.legend()
plt.title("Tracé d'Eyring")
plt.show()