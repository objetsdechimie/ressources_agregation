import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

data = pd.read_csv("loi_des_aires.csv")


# fichier avec colonnes T, X, Y, delta_aire  
T = data["T"].values
X = data["X"].values
Y = data["Y"].values
delta_aire = data["deltaaire"].values/2  # incertitude supplémentaire spécifique à chaque point

# param fixes
uX = 0.002  # incertitude sur X
uY = 0.002  # incertitude sur Y

#  Calcul de l'aire 
A = X[:-1]*Y[1:] - X[1:]*Y[:-1]

# Temps moyen 
T_milieu = 0.5 * (T[:-1] + T[1:])

#  Calcul de l’incertitude sur A pour chaque intervalle 

deltaA_segment = 0.5 * (delta_aire[:-1] + delta_aire[1:])

uA = np.sqrt(
    (Y[1:] * uX)**2 +
    (X[:-1] * uY)**2 +
    (Y[:-1] * uX)**2 +
    (X[1:] * uY)**2 +
    deltaA_segment**2
)

#  Ajustement d'une droite affine 
def droite(T, a, b):
    return a*T + b

params, cov = curve_fit(droite, T_milieu, A, sigma=uA, absolute_sigma=True)
a, b = params
ua, ub = np.sqrt(np.diag(cov))


T_fit = np.linspace(min(T_milieu), max(T_milieu), 200)
A_fit = droite(T_fit, a, b)

# tracé
plt.errorbar(T_milieu, A, yerr=uA, fmt='o', capsize=3, label="Données expérimentales")
plt.plot(T_fit, A_fit, 'r-')
plt.xlabel("T (ms)")
plt.ylabel("Aire A (m$^2$)")
plt.title("Évolution de l'aire en fonction du temps ")
plt.legend()
plt.grid(True)
plt.show()

