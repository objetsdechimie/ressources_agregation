# -*- coding: utf-8 -*-
"""
Created on Thu Mar 26 09:27:33 2026

@author: maxim
"""

import numpy as np
import matplotlib.pyplot as plt

# Donnees initiales
v0 = np.array([0.42, 0.97, 1.71, 2.49, 3.53, 3.72])  # Âµmol/L/min
C0 = np.array([2, 5, 10, 20, 50, 100])  # mmol/L

# Inverses pour Lineweaver-Burk
inv_C0 = 1 / C0
inv_v0 = 1 / v0

# Ajustement lineaire
coeffs = np.polyfit(inv_C0, inv_v0, 1)  # coeffs[0] = a, coeffs[1] = b
a, b = coeffs

# Calcul de Vmax et Km
Vmax = 1 / b
Km = a * Vmax

# Affichage des resultats
print(f"Vmax = {Vmax:.2f} Âµmol/L/min")
print(f"Km = {Km:.2f} mmol/L")

# Trace du graphe avec droite ajustee
x_fit = np.linspace(0, max(inv_C0) * 1.1, 100)
y_fit = a * x_fit + b

plt.figure()
plt.plot(inv_C0, inv_v0, 'ro', label='Donnees experimentales')
plt.plot(x_fit, y_fit, 'b--', label='Ajustement lineaire')
plt.title("Droite de Lineweaver-Burk")
plt.xlabel('1 / [S]$_0$ (L/mmol)')
plt.ylabel('1 / V$_0$ (minÂ·L/Âµmol)')
plt.legend()
plt.grid(True)
plt.show()

# Evolution des concentrations

from scipy.integrate import odeint
import numpy as np
import matplotlib.pyplot as plt

# Parametres cinetiques
k1 = 10     # (L/Âµmol/min)
k_minus1 = 5  # (/min)
k2 = 0.2     # (/min)

# Conditions initiales
S0 = 10     # Âµmol/L
E0 = 1       # Âµmol/L
ES0 = 0      # Âµmol/L
P0 = 0       # Âµmol/L

# Systeme d'equations differentielles
def d_system(y, t, k1, k_minus1, k2):
    S, E, ES, P = y
    dSdt = -k1 * E * S + k_minus1 * ES
    dEdt = -k1 * E * S + (k_minus1 + k2) * ES
    dESdt = k1 * E * S - (k_minus1 + k2) * ES
    dPdt = k2 * ES
    return [dSdt, dEdt, dESdt, dPdt]

# Temps de simulation
t = np.linspace(0, 75, 200)  # minutes

# Resolution
y0 = [S0, E0, ES0, P0]
sol = odeint(d_system, y0, t, args=(k1, k_minus1, k2))
S, E, ES, P = sol.T

# Trace
plt.figure(figsize=(10, 6))
plt.plot(t, S, label='[S] Substrat', color='orange')
plt.plot(t, E, label='[E] Enzyme libre', color='blue')
plt.plot(t, ES, label='[ES] Complexe enzyme-substrat', color='purple')
plt.plot(t, P, label='[P] Produit', color='green')
plt.title("evolution des concentrations : modele enzymatique complet")
plt.xlabel("Temps (min)")
plt.ylabel("Concentration (Âµmol/L)")
plt.legend()
plt.grid(True)
plt.show()