# -*- coding: utf-8 -*-
"""
Created on Tue Mar 17 10:27:50 2026

@author: ens
"""

import matplotlib
import matplotlib.pyplot as plt
import scipy 
import scipy.constants as cst
import numpy as np 

# =============================================================================
# Importation des variables
# =============================================================================

R = cst.R # constante des gaz parfaits, R = 8,314 J/K/mol
F = 96485 # constante de Faraday, en C/mol
T = 293 # température ambiante, en K

p = -(R*T*np.log(10))/F

n = np.array([0,1,2,3,4,5,6,7,8]) # à remplir
E = np.array([-98,-105,-108,-111,-113,-116,-116,-119,-121])*1e-3 # à remplir

y = 10**(E/p)

# =============================================================================
# Régression linéaire
# =============================================================================

coeffs = np.polyfit(n, y, deg=1) 
a, b = coeffs
x = np.linspace(-10,15)

yth = a*x + b  

# =============================================================================
# Tracer des données
# =============================================================================

plt.plot(n, y, 'x', color = 'crimson')
plt.plot(x, yth, '--', color='cornflowerblue' )
plt.plot(x, np.zeros(len(x)), ':', color='grey')
plt.title('Méthode des ajouts dosés')
plt.xlabel("Nombre d'ajouts")
plt.ylabel(r'$10^{E/p}$')
plt.text(-10,250, r'$a = {:.2f}, b = {:.2f}$'.format(a,b) 
         + '\n' r'$n_0 = {:.2f}$'.format(-b/a))
plt.show()

print('n_0 =', -b/a)