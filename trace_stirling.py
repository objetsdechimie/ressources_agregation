# -*- coding: utf-8 -*-
"""
Created on Fri Nov 28 19:09:56 2025

@author: maxim
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

df = pd.read_csv("scope_0_sans_frottements.csv")
df.columns = ["time_s", "volume_V", "pressure_V"]

temps = list(df['time_s'])[1:655]
temps = [float(v.replace(',', '.')) for v in temps]
volume = list(df['volume_V'])[1:655]
volume = [(32+float(v.replace(',', '.'))*12/4.096)*1e-6 for v in volume]
pression = list(df['pressure_V'])[1:655]
pression = [float(v.replace(',', '.'))/20e-6+1e5 for v in pression]

# Conversion en arrays numpy
V = np.array(volume)
P = np.array(pression)

# Calcul du travail avec la formule du lacet (Shoelace formula)
# W = ∮P dV = 1/2 * Σ(P_i + P_{i+1}) * (V_{i+1} - V_i)
travail = 0
for i in range(len(V)-1):
    travail += (P[i] + P[i+1]) * (V[i+1] - V[i]) / 2

print(f"Travail du moteur : {travail:.3f} J")
print(f"Travail du moteur : {travail*1000:.3f} mJ")
# Durée d'un cycle
duree_cycle = temps[652] - temps[1]  # en secondes

# Puissance moyenne
puissance = travail / duree_cycle  # en Watts

print(f"Travail par cycle : {travail:.3f} J")
print(f"Durée du cycle : {duree_cycle:.3f} s")
print(f"Puissance moyenne : {puissance:.3f} W")
# Visualisation
plt.figure(figsize=(10, 6))
plt.plot(V, P, 'b-', linewidth=2)
plt.fill(V, P, alpha=0.3, label=f'Travail = {travail:.3f} J | Puissance = {puissance:.3f} W \n  Temps de cycle = {duree_cycle:.3f} s')
plt.xlabel('Volume (m³)')
plt.ylabel('Pression (Pa)')
plt.title('Diagramme de Clapeyron')
plt.grid(True)
plt.legend()
plt.show()