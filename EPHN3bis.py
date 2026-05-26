#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 25 14:53:49 2026

@author: Agregatif_Lyon assisted by Gemini
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Traitement de données Électrochimie : Argent-Ammoniac
Régression de Monte-Carlo, Khi-deux et Ellipses d'incertitude
"""

import numpy as np
import matplotlib.pyplot as plt 
import pandas as pd
from matplotlib.patches import Ellipse

# --- 1. PARAMÈTRES ET INCERTITUDES TYPES (LOI RECTANGULAIRE/TRIANGULAIRE) ---
u_V = 0.04e-3 / np.sqrt(6)   # Burette graduée (mL -> L)
u_V0 = 0.05e-3 / np.sqrt(6)  # Pipette jaugée (mL -> L)
u_E = 0.0001                 # Voltmètre (0.1 mV)
u_pH = 0.01                  # pH-mètre (centièmes)

# --- 2. LECTURE DES DONNÉES ---
try:
    # 1. Lecture ultra-flexible du séparateur (virgule, point-virgule, tabulation...)
    df = pd.read_csv('E_PNH3_bibi.csv', sep=None, engine='python', dtype=str)

    # 2. Nettoyage généraliste (on gère le cas où les nombres sont vus comme du texte)
    # On remplace les virgules éventuelles par des points et on convertit en float
    def clean_column(series):
        return pd.to_numeric(series.str.replace(',', '.'), errors='coerce')

    V_raw = clean_column(df.iloc[:, 0])
    pH_raw = clean_column(df.iloc[:, 1])
    E_raw  = clean_column(df.iloc[:, 2])

    # 3. Suppression des lignes vides ou corrompues (NaN)
    mask_valid = V_raw.notna() & pH_raw.notna() & E_raw.notna()
    
    V = V_raw[mask_valid].values * 1e-3  # mL en L
    pH = pH_raw[mask_valid].values
    E = E_raw[mask_valid].values

    if len(V) == 0:
        raise ValueError("Aucune donnée numérique valide n'a pu être extraite.")

    print(f"Lecture réussie : {len(V)} points importés.")

except Exception as e:
    print(f"❌ Erreur critique de lecture : {e}")
    print("Conseil : Enregistrez votre fichier en 'CSV (séparateur point-virgule)'.")
    raise

# --- 3. CONSTANTES PHYSIQUES ---
V0, n0Ag, n0NH4 = 50E-3, 5E-4, 2E-3
Ka, E0Ag = 10**(-9.25), 0.8
R, T, F = 8.314, 298, 96500
Beta1 = 10**3.4
nernst = (R * T * np.log(10)) / F # ~ 0.0592 V

# --- 4. FONCTION DE CALCUL ET PROPAGATION ---

def compute_pNH3(V_in, pH_in, E_in, E0_ref):
    """
    Calcule le pNH3 (logarithme de la concentration en ammoniac libre)
    en propageant les mesures expérimentales à travers les équations de Nernst 
    et de conservation de la matière.
    """
    # Force les entrées en tableaux NumPy pour permettre le calcul vectoriel (sur toute la colonne)
    V_in, pH_in, E_in = np.atleast_1d(V_in), np.atleast_1d(pH_in), np.atleast_1d(E_in)
    
    # Calcul du volume total à chaque ajout de la burette (V0 = volume initial dans le bécher)
    V_tot = V_in + V0
    
    # Calcul du potentiel corrigé (Ecorr) : 
    # On normalise le potentiel mesuré par rapport au premier point (E0_ref) 
    # et on intègre le terme de Nernst lié à la concentration initiale en Ag+
    Ec = E_in + E0Ag + (R*T/F) * np.log(n0Ag/V0) - E0_ref
    
    # Terme exponentiel issu de la loi de Nernst : correspond au rapport [Ag+]/[Ag_initial]
    # Inversion de la formule : E = E0 + (RT/F)*ln([Ag+]) => [Ag+] = exp((E-E0)*F/RT)
    t_exp = np.exp((Ec - E0Ag) * F / (R * T))
    
    # NUMÉRATEUR : Bilan de matière sur l'ammoniac total introduit diminué de ce qui 
    # est consommé par la formation du complexe prédominant [Ag(NH3)2]+
    # On multiplie par 2 car chaque ion Ag+ complexe deux molécules de NH3
    num = (n0NH4/V_tot) - 2*(n0Ag/V_tot) + 2*t_exp
    
    # DÉNOMINATEUR : Prise en compte de l'équilibre acido-basique du solvant (NH4+/NH3)
    # et de la stabilité du premier complexe [Ag(NH3)]+ (via Beta1)
    # Terme 10**(-pH)/Ka représente la proportion de NH4+ par rapport à NH3
    den = 1 + (10**(-pH_in)/Ka) - Beta1*t_exp
    
    # Retourne le pNH3 = -log10([NH3])
    # On utilise np.abs pour éviter les erreurs numériques (log de nombre négatif) dues au bruit
    return -np.log10(np.abs(num/den))

# --- CALCUL DES GRANDEURS D'INTÉRÊT ---

# Calcul de Ecorr : potentiel corrigé
Ecorr = E + E0Ag + (R*T/F) * np.log(n0Ag/V0) - E[0]

# Calcul du vecteur pNH3
pNH3 = compute_pNH3(V, pH, E, E[0])

# --- PROPAGATION DES INCERTITUDES (MÉTHODE GUM / LPI) ---

# 1. Incertitude sur Ecorr
# Comme Ecorr = E - E[0] + constante, l'incertitude est sqrt(u_E**2 + u_E**2)
u_Ecorr = np.full_like(Ecorr, np.sqrt(2 * u_E**2))

# 2. Incertitude sur pNH3 par dérivées numériques (Sensibilités)
eps = 1e-6 # Petit incrément pour simuler la dérivée

# Calcul des sensibilités (pente de la fonction pNH3 par rapport à chaque variable)
dp_dV  = (compute_pNH3(V + eps, pH, E, E[0]) - pNH3) / eps
dp_dpH = (compute_pNH3(V, pH + eps, E, E[0]) - pNH3) / eps
dp_dE  = (compute_pNH3(V, pH, E + eps, E[0]) - pNH3) / eps

# Combinaison quadratique des incertitudes (Loi de Propagation des Incertitudes)
# On pondère chaque incertitude-type par sa sensibilité au carré
u_pNH3 = np.sqrt((dp_dV * u_V)**2 + (dp_dpH * u_pH)**2 + (dp_dE * u_E)**2)

# Désormais, tu as deux vecteurs u_Ecorr et u_pNH3 prêts pour tes ellipses et tes régressions.

# --- 5. RÉGRESSION MONTE-CARLO & KHI-DEUX ---

def regression_monte_carlo(x_global, y_global, ux_global, uy_global, x_min, x_max):
    N_SIMULATIONS = 10000   
    T_STUDENT = 2.571  # Pour n=7 points env. Ajuster si besoin.

    mask = (x_global >= x_min) & (x_global <= x_max) & np.isfinite(x_global) & np.isfinite(y_global)
    X, sX, Y, sY = x_global[mask], ux_global[mask], y_global[mask], uy_global[mask]

    if len(X) < 2: return None

    # Régression pondérée initiale
    a_init, b_init = np.polyfit(X, Y, 1)
    p_weights = 1.0 / np.sqrt(sY**2 + (a_init * sX)**2)
    a_khi2, b_khi2 = np.polyfit(X, Y, 1, w=p_weights)

    # Monte-Carlo
    rng = np.random.default_rng()
    X_rd = rng.normal(X, sX, (N_SIMULATIONS, len(X)))
    Y_rd = rng.normal(Y, sY, (N_SIMULATIONS, len(Y)))
    
    a_rd, b_rd = [], []
    for i in range(N_SIMULATIONS):
        p_i = 1.0 / np.sqrt(sY**2 + (a_khi2 * sX)**2)
        coeffs = np.polyfit(X_rd[i], Y_rd[i], 1, w=p_i)
        a_rd.append(coeffs[0])
        b_rd.append(coeffs[1])

    s_a, s_b = np.std(a_rd), np.std(b_rd)
    return {"a": a_khi2, "U_a": T_STUDENT * s_a, "b": b_khi2, "U_b": T_STUDENT * s_b, "X": X}

# --- 6. ELLIPSES ET AFFICHAGE ---

def draw_error_ellipse(ax, x, y, sx, sy, color='lightgray'):
    ellipse = Ellipse((x, y), width=2*sx, height=2*sy, edgecolor=color, facecolor='none', alpha=0.5)
    ax.add_patch(ellipse)

fig, ax = plt.subplots(figsize=(12, 8))

# Ellipses et points
for i in range(len(pNH3)):
    if np.isfinite(pNH3[i]) and np.isfinite(Ecorr[i]):
        draw_error_ellipse(ax, pNH3[i], Ecorr[i], u_pNH3[i], u_Ecorr[i])

ax.errorbar(pNH3, Ecorr, xerr=u_pNH3, yerr=u_Ecorr, fmt='k.', ecolor='gray', alpha=0.6, label='Données exp.')

# Régressions
# --- CALCUL ZONE 1 ---
res1 = regression_monte_carlo(pNH3, Ecorr, u_pNH3, u_Ecorr, 2, 3)

# Tracé Zone 1
x_plot1 = np.linspace(2, 3, 10)
ax.plot(x_plot1, res1["a"] * x_plot1 + res1["b"], color='red', lw=2, 
        label=f'Zone 1 : pente={res1["a"]:.3f} ± {res1["U_a"]:.3f}')

# Affichage stoechiométrie
print(f"Stoechiométrie p = {abs(res1['a'])/nernst:.2f} ± {res1['U_a']/nernst:.2f}")


# --- CALCUL ZONE 2 ---
res2 = regression_monte_carlo(pNH3, Ecorr, u_pNH3, u_Ecorr, 4.5, 6)

# Tracé Zone 2
x_plot2 = np.linspace(4.5, 6, 10)
ax.plot(x_plot2, res2["a"] * x_plot2 + res2["b"], color='green', lw=2, 
        label=f'Zone 2 : pente={res2["a"]:.3f} ± {res2["U_a"]:.3f}')
            
            

ax.set_xlabel('$pNH_3$')
ax.set_ylabel('$E_{corr}$ (V)')
ax.set_title('Diagramme Potentiel - $pNH_3$ (Monte-Carlo & Ellipses)')
ax.legend()
ax.grid(True, linestyle=':', alpha=0.6)
plt.show()