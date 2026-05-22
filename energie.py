# -*- coding: utf-8 -*-
"""
Created on Thu Mar  5 09:07:59 2026

@author: maxim
"""


import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button # outils pour faire des curseurs

alpha=-8.1 # Energie arbitraire de l'orbitale consideree
beta= -1 # Valeur de l'integrale de recouvrement entre deux orbitales voisines (arbitraire)
n_max=1000
a=1e-10
kb = 1.35e-23
T=298

def E(N,n) :
    # Renvoie l'energie du n-ieme niveau d'energie d'une chaÃ®ne de N atomes
    return alpha+2*beta*np.cos(n*np.pi/(N+1))

Emin = alpha + 2*beta

Emax = alpha - 2*beta

fig,ax=plt.subplots()
A=[0]*n_max
N_init=1 # Nombre d'atomes initialement
s=0
for i in range(0,n_max):
    if i < N_init:
        Enew = E(N_init, i+1)
        if i == 0:
            A[i], = plt.plot([-i,i],[Enew,Enew],color='red',label='libre')
        else:
            A[i], = plt.plot([-i,i],[Enew,Enew],color='red')
    else:
        if i == N_init:
            A[i], = plt.plot([-i,i],[Enew+10,Enew+10],color='black',label='occupé')
        else:
            A[i], = plt.plot([-i,i],[Enew+10,Enew+10],color='black')
    
        
ax.set_ylim(-10.4,-5.8)
ax.set_ylabel('Energie (eV)',fontsize=15)
ax.set_title('Energie des niveaux electroniques (methode de Huckel)',fontsize=20)
plt.xticks(color='white')
plt.yticks(fontsize=13)
plt.legend(loc='upper right', )

# Ajuste le graphique pour faire de la place pour le curseur
plt.subplots_adjust(bottom=0.25)

axN = plt.axes([0.25, 0.1, 0.65, 0.03])
N_slider = Slider(
    ax=axN,
    label='Nombres d\'atomes',
    valmin=1,
    valmax=n_max,
    valstep=1,
    valinit=1)

def fermi(E,ef,T,kb):

    return 1/(1+np.exp((E-ef)/(kb*T)))

# Fonction appelee des que la valeur de theta change
def update(val):
    ef=val//2 +1
    for i in range(0,n_max) :#Trace des N_init niveaux d'energies
        if ef < i :
            if i<N_slider.val :
                Enew=E(N_slider.val,i+1)
                A[i].set_ydata([Enew,Enew])
                A[i].set_color('red')
            if i>=N_slider.val:
                A[i].set_ydata([0,0])
                A[i].set_color('white')
        else : 
            if i<N_slider.val :
                Enew=E(N_slider.val,i+1)
                A[i].set_ydata([Enew,Enew])
                A[i].set_color('black')
            if i>=N_slider.val:
                A[i].set_ydata([0,0])
                A[i].set_color('white')

N_slider.on_changed(update)
