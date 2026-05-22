# -*- coding: utf-8 -*-
"""
Created on Fri May 15 13:58:31 2026

@author: maxim
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider
from scipy.optimize import fsolve

# --- Configuration de la réaction ---
ORDER = 2
K = 1
C0 = 1
X_TARGET = 0.9

def inv_r(X):
    return 1 / (K * (C0 * (1 - X))**ORDER)

def solve_intermediate_conversions(n_stages, x_final):
    """Calcule les X_i pour n réacteurs de volumes égaux en série."""
    if n_stages == 1:
        return [0, x_final]
    
    # On cherche tau_par_etage tel que la conversion finale soit x_final
    def equations(vars):
        # vars sont [tau, X1, X2, ..., X_{n-1}]
        tau = vars[0]
        conversions = list(vars[1:]) + [x_final]
        res = []
        prev_X = 0
        for current_X in conversions:
            # Bilan : tau = (X_i - X_{i-1}) / r(X_i)
            res.append(tau - (current_X - prev_X) * inv_r(current_X))
            prev_X = current_X
        return res

    # Initialisation : répartition linéaire de X et un tau estimé
    initial_guess = [x_final/n_stages] + list(np.linspace(0, x_final, n_stages+1)[1:-1])
    sol = fsolve(equations, initial_guess)
    
    tau_total = sol[0] * n_stages
    all_X = [0] + list(sol[1:]) + [x_final]
    return all_X

# --- Création du Graphique ---
fig, ax = plt.subplots(figsize=(10, 7))
plt.subplots_adjust(bottom=0.25) # Place pour le curseur

X_plot = np.linspace(0, 0.95, 200)
Y_plot = inv_r(X_plot)
line, = ax.plot(X_plot, Y_plot, 'k-', lw=2, label='1/(-r)')

# Zone RP (Piston) - Fixe
X_rp = np.linspace(0, X_TARGET, 100)
ax.fill_between(X_rp, inv_r(X_rp), color='green', alpha=0.2, label='Volume RP')

# Conteneur pour les rectangles des RPAC
rects = []

def update_plot(val):
    n = int(val)
    # Nettoyer les anciens rectangles
    for r in rects:
        r.remove()
    rects.clear()
    
    # Calculer les nouvelles étapes
    conversions = solve_intermediate_conversions(n, X_TARGET)
    
    # Dessiner les nouveaux rectangles
    total_tau = 0
    for i in range(len(conversions)-1):
        x_start = conversions[i]
        x_end = conversions[i+1]
        h = inv_r(x_end)
        width = x_end - x_start
        total_tau += width * h
        
        rect = ax.bar(x_start, h, width=width, align='edge', 
                      alpha=0.4, color='royalblue', edgecolor='blue', lw=1)
        rects.append(rect[0])
    
    ax.set_title(f'Série de {n} RPAC | Tau Total : {total_tau:.2f} (RP : {4.5:.2f} env.)')
    fig.canvas.draw_idle()

# --- Curseur ---
ax_slider = plt.axes([0.2, 0.1, 0.6, 0.03])
slider_n = Slider(ax_slider, 'Nb de RPAC', 1, 20, valinit=1, valstep=1)
slider_n.on_changed(update_plot)

# Initialisation
update_plot(1)

ax.set_xlim(0, 0.95)
ax.set_ylim(0, 110)
ax.set_xlabel('Conversion (X)')
ax.set_ylabel('1 / (-r)')
ax.legend(loc='upper left')
ax.grid(True, alpha=0.3)

plt.show()