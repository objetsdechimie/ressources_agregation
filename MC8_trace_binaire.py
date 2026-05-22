#!/bin/env python3

# from math import *
import matplotlib.pyplot as plt
import numpy as np
import scipy.optimize

"""
Date : 2019
Auteur : Joachim Galiana, Manon Leconte


Ce programme vise à tracer les branches d'un diagramme binaire solide/liquide possédant un eutectique.
Les deux branches sont calculées et tracées et leur intersection (eutectique) est calcuée et pointée de même.
"""

#DÉFINITION DES GRANDEURS et CONSTANTES PHYSIQUES
Tfus_d=352.4#en K
Tfus_p=372.4#en K
R=8.314#en J/K/mol
Delta_Fus_d=21e3
Delta_Fus_p=16.46e3
M_d=134.2
M_p=178.2
mtot=2.0

#DÉFINITION DES VARIABLES
x=np.linspace(0.01,0.99,99)
y=1-x

#AIDE AU TRACÉ EXPÉRIMENTAL
x_exp=np.array([0,0.1,0.2,0.3,0.4,0.54,0.6,0.7,0.8,0.9,1])
print(x_exp)
m_p=mtot/(1 + (M_d * x_exp)/(M_p * (1 - x_exp)))
m_d=2-m_p
print(m_d)
print(m_p)


#DÉFINITION DES FONCTIONS
def Td(fractions):
    return((Tfus_d*Delta_Fus_d)/(Delta_Fus_d-R*Tfus_d*np.log(fractions)))
def Tp(fractions):
    return((Tfus_p*Delta_Fus_p)/(Delta_Fus_p-R*Tfus_p*np.log(fractions)))
def Soustraction(fractions):
    return (Tfus_d*Delta_Fus_d)/(Delta_Fus_d-R*Tfus_d*np.log(fractions)) - (Tfus_p*Delta_Fus_p)/(Delta_Fus_p-R*Tfus_p*np.log(1-fractions))
# T1=np.array([(Td*Delta_Fus_d)/(Delta_Fus_d-R*Td*log(a)) for a in x])
# T2=np.array([(Tp*Delta_Fus_p)/(Delta_Fus_p-R*Tp*log(b)) for b in y])

# Tp=Tp(x)
# Td=Td(y)

"""
L'abcisse utilisée est x, fraction molaire en Durène : c'est toujours par rapport au Durène qu'on réfléchit à partir d'ici.
"""

E=scipy.optimize.fsolve(Soustraction,0)
print(E)
ante_x=np.array([a for a in x if a < E])
post_x=np.array([a for a in x if a > E])

plt.ylim(285,380)
# plt.plot(x,Td(x),"red")
plt.plot(ante_x,Td(ante_x),"red",linestyle='--')
plt.plot(post_x,Td(post_x),"red")
plt.plot(ante_x,Tp(1-ante_x),"blue")
plt.plot(post_x,Tp(1-post_x),"blue",linestyle='--')
# plt.plot(x,Tp(y))
plt.plot(E,Td(E),"+")
plt.hlines(y=Td(E),xmin=0,xmax=1)
# plt.hlines(x=E[0],ymin=0,ymax=Td(E)[0])
# plt.plot(x,Soustraction(x))
plt.xlabel("Fraction molaire en durène, x")
plt.ylabel("Température, T (K)")
plt.title("Diagramme binaire solide/liquide durène-phénanthrène\n Eutectique à x = 0,54")
plt.grid()
plt.show()
