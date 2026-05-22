# -*- coding: utf-8 -*-
"""
Created on Thu Mar 26 09:35:18 2026

@author: maxim
"""

import numpy as np
import matplotlib.pyplot as plt

N = 100000 # Nombre de simulations pour la methode de Monte-Carlo

# Definition du volume de la solution titree
V0 = 25 # Volume initial de la solution titree (en mL)
EMT_V0 = 0.03 # Erreur maximale toleree sur la pipette (en mL)

# Definition du volume a lâ€™equivalence
Veq = 19.40 # Volume a lâ€™equivalence (en mL)
EMT_burette = 0.03 # Erreur maximale toleree sur la burette (en mL)
graduation = 0.025 # Demi-graduation de la burette (en mL)
goutte_determination = 0.05 # Determination de l'equivalence a la goutte pres (en mL)

# Definition de la concentration de la solution titrante preparee a partir d'une pesee
m = 0.0238 # Masse de solide du titrant (en g)
EMT_m = 0.0001 # Erreur maximale toleree sur la masse de solide (en g)
Vfiole = 100e-3 # Volume de la fiole (en L)
EMT_Vfiole = 0.1e-3 # Erreur maximale toleree sur le volume de la fiole (en L)
M = 158.11 # Masse molaire du titrant (en g/mol)

# Methode de Monte-Carlo
V0_rd = np.random.triangular(V0-EMT_V0, V0, V0+EMT_V0, N)
Veq_rd = Veq + np.random.triangular(-EMT_burette, 0, EMT_burette, N) + np.random.uniform(-graduation, graduation, N) + np.random.uniform(-graduation, graduation, N) + np.random.uniform(-goutte_determination, goutte_determination, N)
m_rd = m + np.random.uniform(-EMT_m, EMT_m, N) + np.random.uniform(-EMT_m, EMT_m, N) 
Vfiole_rd = np.random.triangular(Vfiole-EMT_Vfiole, Vfiole, Vfiole+EMT_Vfiole, N)
C = (32*m_rd*Veq_rd)/(4*V0_rd*M*Vfiole_rd)*1e3 # Calcul pour determiner la concentration massique en dioxygene en mg/L
C_calc = (32*m*Veq)/(4*V0*M*Vfiole)*1e3
s_c = np.std(C)

# Affichage du resultat
print('c(O2) =', format(C_calc,'.2E'),'\u00B1', format(2*s_c,'.1E'), 'mg/L a 95 % de confiance')

# Histogramme des differentes valeurs de C
plt.figure("Tirage")
plt.hist(C, bins = 200, color = 'grey', alpha = 0.7)
plt.xlabel("c(O$_2$) (en mg/L)", fontsize = 20)
plt.ylabel("Occurrences", fontsize = 20)
plt.xticks(fontsize = 20)
plt.yticks(fontsize = 20)
plt.show()

# Diagramme en barres de la contribution relative des incertitudes
C_V0 = (32*m*Veq)/(4*V0_rd*M*Vfiole)*1e3 # Creation d'un tableau de concentration avec uniquement V0 qui varie
C_Veq = (32*m*Veq_rd)/(4*V0*M*Vfiole)*1e3 # Creation d'un tableau de concentration avec uniquement Veq qui varie
C_m = (32*m_rd*Veq)/(4*V0*M*Vfiole)*1e3 # Creation d'un tableau de concentration avec uniquement m qui varie
C_Vfiole = (32*m*Veq)/(4*V0*M*Vfiole_rd)*1e3 # Creation d'un tableau de concentration avec uniquement Vfiole qui varie
simul = [C_V0, C_Veq, C_m, C_Vfiole]

valeur =[]
prop_tot = np.var(C_V0) + np.var(C_Veq) + np.var(C_m) + np.var(C_Vfiole) # Variance totale
for i in range(len(simul)):
     valeur.append(np.var(simul[i])/prop_tot*100) # Calcul de la contribution relative a la variance
 
plt.figure("Proportion")
plt.bar([1,2,3,4], height = valeur, tick_label = ('$V_0$', '$V_{eq}$', '$m$', '$V_{fiole}$'), color = 'grey', alpha = 0.7)
plt.xlabel("Grandeur mesuree", fontsize = 20)
plt.ylabel("Contribution relative de chaque incertitude (%)", fontsize = 20)
for i in range(len(valeur)):
        plt.text(i+1, valeur[i]+1, round(valeur[i],1), ha = 'center', fontsize = 20)
plt.xticks(fontsize = 20)
plt.yticks(fontsize = 20)
plt.show()