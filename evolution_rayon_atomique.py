# -*- coding: utf-8 -*-
"""
Created on Sat Mar 21 11:57:59 2026

@author: maxim
"""

import matplotlib.pyplot as plt

# ==========================================================
# 1. MOTEUR DE CALCUL DE SLATER
# ==========================================================

def get_slater_radius(Z, config, target_group):
    """Calcule le rayon atomique théorique selon les règles de Slater."""
    # Groupes de Slater dans l'ordre : (nom, n, type)
    groups = [
        ('1s', 1, 's'), ('2s2p', 2, 'sp'), ('3s3p', 3, 'sp'), 
        ('3d', 3, 'd'), ('4s4p', 4, 'sp'), ('4d', 4, 'd'), 
        ('4f', 4, 'f'), ('5s5p', 5, 'sp'), ('5d', 5, 'd'), 
        ('6s6p', 6, 'sp')
    ]
    
    # Identification du groupe cible
    t_n, t_type = next((n, gtype) for name, n, gtype in groups if name == target_group)
    
    S = 0.0
    for g_name, g_n, g_type in groups:
        count = config.get(g_name, 0)
        if count == 0: continue
        
        # On retire l'électron sur lequel on calcule la force
        adj_count = count - 1 if g_name == target_group else count
        
        if g_name == target_group:
            S += adj_count * (0.30 if g_name == '1s' else 0.35)
        elif g_n == t_n: # Même n mais à gauche (ex: 3d vs 4s)
            S += adj_count * 0.35
        elif t_type == 'sp':
            if g_n == t_n - 1: S += adj_count * 0.85
            elif g_n < t_n - 1: S += adj_count * 1.0
        else: # Cible est d ou f
            S += adj_count * 1.0
            
    Zeff = Z - S
    # Correction de n selon Slater (n effectif)
    n_star = {1:1, 2:2, 3:3, 4:3.7, 5:4, 6:4.2}[t_n]
    return 0.529 * (n_star**2 / Zeff)

# ==========================================================
# 2. GÉNÉRATION DES DONNÉES (LIGNES 1 À 4)
# ==========================================================

symbols = [
    "H","He","Li","Be","B","C","N","O","F","Ne","Na","Mg","Al","Si","P","S","Cl","Ar",
    "K","Ca","Sc","Ti","V","Cr","Mn","Fe","Co","Ni","Cu","Zn","Ga","Ge","As","Se","Br","Kr"
]

data_main = []
for i, sym in enumerate(symbols):
    Z = i + 1
    # Détermination simplifiée de la config pour Slater
    conf = {'1s': min(Z, 2)}
    if Z > 2:  conf['2s2p'] = min(Z-2, 8)
    if Z > 10: conf['3s3p'] = min(Z-10, 8)
    if Z > 18: 
        # Remplissage 4s puis 3d (simplifié pour Slater)
        if Z <= 20: 
            conf['4s4p'] = Z-18
        else:
            conf['3d'] = min(Z-20, 10)
            conf['4s4p'] = max(2, Z-28) # Electrons p après le bloc d
            
    # Déterminer la couche de valence
    v_layer = '1s' if Z <= 2 else '2s2p' if Z <= 10 else '3s3p' if Z <= 18 else '4s4p'
    data_main.append({'Z': Z, 'sym': sym, 'r': get_slater_radius(Z, conf, v_layer)})

# ==========================================================
# 3. GÉNÉRATION DES DONNÉES (LANTHANIDES)
# ==========================================================

lanth_syms = ["La","Ce","Pr","Nd","Pm","Sm","Eu","Gd","Tb","Dy","Ho","Er","Tm","Yb","Lu"]
data_lanth = []
for i, sym in enumerate(lanth_syms):
    Z = 57 + i
    # Config [Xe] 4f^N 6s^2 -> Slater : (1-5p) 4f^N 6s^2
    # Le groupe [Xe] complet = 54 e-
    conf = {'1s':2, '2s2p':8, '3s3p':8, '3d':10, '4s4p':8, '4d':10, '5s5p':8, '4f': i, '6s6p':2}
    data_lanth.append({'Z': Z, 'sym': sym, 'r': get_slater_radius(Z, conf, '6s6p')})

# ==========================================================
# 4. TRACÉ DES GRAPHIQUES
# ==========================================================

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(18, 7))

# --- Graphe 1 : Périodes 1-4 ---
colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728']
for p in range(1, 5):
    z_range = {1:(1,2), 2:(3,10), 3:(11,18), 4:(19,36)}[p]
    subset = [d for d in data_main if z_range[0] <= d['Z'] <= z_range[1]]
    zs = [d['Z'] for d in subset]
    rs = [d['r'] for d in subset]
    ax1.plot(zs, rs, '-o', label=f'Période {p}', color=colors[p-1], markersize=4)
    for d in subset:
        ax1.annotate(d['sym'], (d['Z'], d['r']), textcoords="offset points", xytext=(0,5), ha='center', fontsize=8)

ax1.set_title("Rayons des 4 premières lignes (Modèle de Slater)")
ax1.set_xlabel("Numéro Atomique (Z)")
ax1.set_ylabel("Rayon de valence ($Å$)")
ax1.legend()
ax1.grid(True, alpha=0.3)

# --- Graphe 2 : Lanthanides ---
lz = [d['Z'] for d in data_lanth]
lr = [d['r'] for d in data_lanth]
ax2.plot(lz, lr, 'g-s', label="Rayon 6s (Slater)")
for d in data_lanth[::2]: # Un nom sur deux pour la clarté
    ax2.annotate(d['sym'], (d['Z'], d['r']), textcoords="offset points", xytext=(0,5), ha='center', fontsize=9)

ax2.set_title("Stagnation des Lanthanides chez Slater")
ax2.set_xlabel("Numéro Atomique (Z)")
ax2.set_ylabel("Rayon orbitale $6s$ ($Å$)")
ax2.set_ylim(min(lr)-0.2, max(lr)+0.2)
ax2.annotate("Rayon strictement constant !\nCar dZ = dS (Écran 4f = 1.0)", xy=(64, lr[0]), 
             xytext=(58, lr[0]+0.1), arrowprops=dict(arrowstyle="->", color='red'))
ax2.grid(True, alpha=0.3)

plt.tight_layout()
plt.show()