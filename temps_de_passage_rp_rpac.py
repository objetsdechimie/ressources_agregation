# -*- coding: utf-8 -*-
"""
Created on Fri May 15 13:52:09 2026

@author: maxim
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import quad
from scipy.optimize import fsolve

def solve_reactors(order, X_target, n_cstr=1):
    """
    Calcule les temps de passage pour RP, RPAC unique et N RPAC en série.
    On pose k = 1 et C0 = 1 pour normaliser les résultats.
    """
    
    # 1. Temps de passage pour le Réacteur Piston (RP)
    # Formule : tau = intégrale de 0 à X de (1 / (1-X)^n)
    func_rp = lambda x: 1 / (1 - x)**order
    tau_rp, _ = quad(func_rp, 0, X_target)
    
    # 2. Temps de passage pour un RPAC unique
    # Formule : tau = X / (1-X)^n
    tau_rpac_unique = X_target / (1 - X_target)**order
    
    # 3. Temps de passage pour N RPAC en série (volumes égaux)
    # On doit résoudre l'équation d'étage en étage
    def cstr_series_system(tau_total, n_stages, target_X):
        tau_i = tau_total / n_stages
        current_C = 1.0  # C/C0 initial
        for _ in range(n_stages):
            # Equation de bilan : C_in - C_out = tau_i * k * C_out^n
            # On cherche C_out tel que : C_out + tau_i * C_out^n - C_in = 0
            obj = lambda c: c + tau_i * (c**order) - current_C
            current_C = fsolve(obj, current_C)[0]
        return (1 - current_C) - target_X

    # On cherche le tau_total qui permet d'atteindre X_target
    tau_rpac_series = fsolve(cstr_series_system, tau_rp, args=(n_cstr, X_target))[0]
    
    return tau_rp, tau_rpac_unique, tau_rpac_series

# --- Paramètres de simulation ---
conversion_target = 0.90  # 90% de conversion
orders = [0.5, 1, 1.5, 2]
n_tanks = 10 # Nombre de RPAC en série

# --- Calculs et Graphique ---
results_rp = []
results_rpac = []
results_series = []

for n in orders:
    t_rp, t_rpac, t_series = solve_reactors(n, conversion_target, n_tanks)
    results_rp.append(t_rp)
    results_rpac.append(t_rpac)
    results_series.append(t_series)

# Plot
x = np.arange(len(orders))
width = 0.25

fig, ax = plt.subplots(figsize=(10, 6))
ax.bar(x - width, results_rpac, width, label='RPAC Unique', color='#e74c3c')
ax.bar(x, results_series, width, label=f'{n_tanks} RPAC en série', color='#f1c40f')
ax.bar(x + width, results_rp, width, label='RP (Piston)', color='#2ecc71')

ax.set_ylabel('Temps de passage adimensionnel (tau)')
ax.set_xlabel('Ordre de la réaction (n)')
ax.set_title(f'Comparaison des réacteurs pour X = {conversion_target*100}%')
ax.set_xticks(x)
ax.set_xticklabels([f'n={n}' for n in orders])
ax.legend()
ax.grid(axis='y', linestyle='--', alpha=0.7)

plt.tight_layout()
plt.show()