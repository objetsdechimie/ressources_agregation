#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May 15 19:19:31 2026

@author: olivieralbrich-sales
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider
from scipy.optimize import fsolve

# --- 1. Paramètres de l'orbite (Ellipse) ---
a = 5.0   # Demi-grand axe
e = 0.6   # Excentricité (assez élevée pour bien voir la différence de vitesse)
b = a * np.sqrt(1 - e**2)  # Demi-petit axe
c = a * e # Distance centre-foyer

# Le foyer (Soleil) est en (0,0). Le centre de l'ellipse est en (-c, 0)
x_foyer, y_foyer = 0, 0

# --- 2. Résolution de l'équation de Kepler ---
# M = E - e*sin(E). On cherche E à partir de l'anomalie moyenne M.
def kepler_equation(E, M, e):
    return E - e * np.sin(E) - M

def get_position_from_time(t):
    """ Calcule la position (x,y) sur l'ellipse pour un temps t donné """
    # M est proportionnel au temps (période T = 2*pi pour simplifier)
    M = t % (2 * np.pi)
    # Résolution numérique pour trouver l'anomalie excentrique E
    E_guess = M
    E_sol = fsolve(kepler_equation, E_guess, args=(M, e))[0]
    
    # Coordonnées par rapport au centre de l'ellipse
    x_centre = a * np.cos(E_sol)
    y_centre = b * np.sin(E_sol)
    
    # Translaté au foyer (Soleil à l'origine)
    return x_centre - c, y_centre

# --- 3. Initialisation de la figure ---
fig, ax = plt.subplots(figsize=(8, 7))
plt.subplots_adjust(bottom=0.25) # Place pour le curseur
ax.set_aspect('equal')
ax.grid(True, linestyle='--', alpha=0.5)
ax.set_title("Pédagogie : Loi des Aires de Kepler ($A_1 = A_2$)", fontsize=14, fontweight='bold')

# Dessin de l'orbite complète
theta = np.linspace(0, 2*np.pi, 200)
# Équation polaire de l'ellipse par rapport au foyer
p = a * (1 - e**2)
r_orbit = p / (1 + e * np.cos(theta))
x_orbit = r_orbit * np.cos(theta)
y_orbit = r_orbit * np.sin(theta)
ax.plot(x_orbit, y_orbit, color='gray', linestyle='--', label='Orbite elliptique')

# Dessin du Soleil (Foyer)
ax.plot(x_foyer, y_foyer, 'yo', markersize=12, label='Soleil (Foyer)')

# Temps initiaux (positions de départ fixes pour les deux planètes)
t1_dep = 0.0          # Proche du périastre (vitesse rapide)
t2_dep = np.pi - 0.5  # Proche de l'apoastre (vitesse lente)

# Éléments graphiques qui seront mis à jour
patch1 = None
patch2 = None
points_planetes, = ax.plot([], [], 'ro', markersize=8)

# Zone de texte pour afficher les aires
text_aires = ax.text(0.05, 0.05, "", transform=ax.transAxes, fontsize=11,
                     bbox=dict(boxstyle="round", facecolor="wheat", alpha=0.5))

# --- 4. Fonction de mise à jour du graphique ---
def update(delta_t):
    global patch1, patch2
    
    # Nettoyage des anciens coloriages d'aires
    if patch1: patch1.remove()
    if patch2: patch2.remove()
    
    # Échantillonnage du temps pour dessiner les arcs balayés
    t_v1 = np.linspace(t1_dep, t1_dep + delta_t, 100)
    t_v2 = np.linspace(t2_dep, t2_dep + delta_t, 100)
    
    # Calcul des coordonnées des arcs
    x1, y1 = [], []
    for t in t_v1:
        pos = get_position_from_time(t)
        x1.append(pos[0])
        y1.append(pos[1])
        
    x2, y2 = [], []
    for t in t_v2:
        pos = get_position_from_time(t)
        x2.append(pos[0])
        y2.append(pos[1])
        
    # Ajout du foyer pour fermer les polygones (secteurs d'aire)
    x1_poly = [x_foyer] + x1 + [x_foyer]
    y1_poly = [y_foyer] + y1 + [y_foyer]
    x2_poly = [x_foyer] + x2 + [x_foyer]
    y2_poly = [y_foyer] + y2 + [y_foyer]
    
    # Coloriage des aires
    patch1 = ax.fill(x1_poly, y1_poly, color='deepskyblue', alpha=0.6, label='Aire 1 (Périastre)')[0]
    patch2 = ax.fill(x2_poly, y2_poly, color='coral', alpha=0.6, label='Aire 2 (Apoastre)')[0]
    
    # Position actuelle des planètes (au bout du temps delta_t)
    p1_actuel = (x1[-1], y1[-1])
    p2_actuel = (x2[-1], y2[-1])
    points_planetes.set_data([p1_actuel[0], p2_actuel[0]], [p1_actuel[1], p2_actuel[1]])
    
    # Calcul mathématique de l'aire balayée (Loi de Kepler : Aire = (A_totale / T) * delta_t)
    # L'aire totale d'une ellipse est pi * a * b. Notre période T est 2*pi.
    aire_theorique = 0.5 * a * b * delta_t
    
    # Mise à jour du texte
    text_aires.set_text(f"Intervalle de temps \u0394t = {delta_t:.2f} s\n"
                        f"Aire Bleue (proche) = {aire_theorique:.3f} ua\u00b2\n"
                        f"Aire Orange (loin)  = {aire_theorique:.3f} ua\u00b2")
    
    # Gestion de la légende pour éviter les doublons à chaque mouvement de curseur
    handles, labels = ax.get_legend_handles_labels()
    by_label = dict(zip(labels, handles))
    ax.legend(by_label.values(), by_label.keys(), loc='upper right')
    
    fig.canvas.draw_idle()

# --- 5. Création du Curseur (Slider) ---
ax_slider = plt.axes([0.2, 0.1, 0.6, 0.04]) # Position du curseur [gauche, bas, largeur, hauteur]
slider_dt = Slider(
    ax=ax_slider,
    label='Durée (\u0394t) ',
    valmin=0.01,
    valmax=2.0,
    valinit=0.5,
    valfmt='%1.2f s',
    color='mediumseagreen'
)

# Liaison du curseur à la fonction de mise à jour
slider_dt.on_changed(update)

# Limites des axes pour bloquer la vue
ax.set_xlim(-a - c - 0.5, a - c + 0.5)
ax.set_ylim(-b - 0.5, b + 0.5)

# Premier affichage au démarrage
update(0.5)

plt.show()