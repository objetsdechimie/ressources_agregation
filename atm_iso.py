import numpy as np
import matplotlib.pyplot as plt


# Données expérimentales (document 1)


# altitude (m)
z_mes = np.array([
     0,  200,  400,  600,  800,
  1000, 1200, 1400, 1600, 1800,
  2000, 2200, 2400, 2600, 2800,
  3000, 3200, 3400, 3600, 3800,
  4000, 4200, 4400, 4600, 4800,
  5000, 5200, 5400, 5600, 5800,
  6000, 6200, 6400, 6600, 6800,
  7000, 8000, 9000,10000,12000,
 14000,16000,18000
])

# pression expérimentale (kPa)
P_mes_kPa = np.array([
101.33, 98.95, 96.61, 94.32, 92.08,
 89.88, 87.72, 85.60, 83.53, 81.49,
 79.50, 77.55, 75.63, 73.76, 71.92,
 70.12, 68.36, 66.63, 64.94, 63.28,
 61.66, 60.07, 58.52, 57.00, 55.51,
 54.05, 52.62, 51.23, 49.86, 48.52,
 47.22, 45.94, 44.69, 43.47, 42.27,
 41.11, 35.65, 30.80, 26.50, 19.40,
 14.17, 10.53,  7.57
])

# Conversion en Pa
P_mes = P_mes_kPa * 1000


# Constantes physiques
#

g = 9.81              # accélération de la pesanteur (m/s²)
R = 8.314             # constante des gaz parfaits (SI)
M = 0.029             # masse molaire de l’air (kg/mol)

# Température supposée constante (modèle isotherme)
T0 = 288.15           # 15°C en Kelvin


# Méthode d’Euler
#

# Pas de calcul
dz = 100

# altitude maximale
z_max = 18000

# tableau des altitudes
z = np.arange(0, z_max + dz, dz)

# tableau pression modèle
P = np.zeros(len(z))

# condition initiale
P[0] = 101330     # pression au niveau du sol (Pa)

# ------------------------------------------------------------
# Boucle d’Euler
# ------------------------------------------------------------
#
# Formule :
#
# P(i+1) = P(i) + dz * dP/dz
#
# avec :
#
# dP/dz = -(M*g)/(R*T0) * P(i)
#
# ------------------------------------------------------------

for i in range(len(z)-1):

    dPdz = -(M * g) / (R * T0) * P[i]

    P[i+1] = P[i] + dz * dPdz

# Conversion en bar pour les graphiques
P_bar = P / 1e5
P_mes_bar = P_mes / 1e5

#
# Graphe pression : modèle + mesures
#

plt.figure(figsize=(8,5))

plt.plot(z, P_bar,
         'b--',
         label='modèle isotherme T = T0')

plt.plot(z_mes, P_mes_bar,
         'ro',
         label='mesures')

plt.xlabel('z (m)')
plt.ylabel('p (bar)')
plt.title('Pression atmosphérique')
plt.grid()

plt.legend()


#  Calcul de l’écart relatif


P_modele_mes = np.interp(z_mes, z, P)


ecart = (P_modele_mes - P_mes) / P_mes


# Graphe de l’écart relatif


plt.figure(figsize=(8,5))

plt.plot(z_mes, ecart, 'b.')

plt.xlabel('z (m)')
plt.ylabel('(p_modèle - p_mesures)/p_mesures')

plt.title('Écart relatif du modèle isotherme')

plt.grid()

plt.show()