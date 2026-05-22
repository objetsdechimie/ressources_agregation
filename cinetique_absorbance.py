# -*- coding: utf-8 -*-
"""
Created on Wed Apr  1 19:32:10 2026

@author: maxim
"""


import matplotlib.pyplot as plt
import numpy as np
from matplotlib.widgets import Slider
import scipy.optimize as spo
import matplotlib
matplotlib.rc('xtick', labelsize=24) 
matplotlib.rc('ytick', labelsize=24) 
matplotlib.rcParams.update({'font.size': 22})

##Valeurs
"""Les valeurs sont issues de Gruber"""

t = np.arange(0,381,20)+14      #Premier set de donnees
A = np.array([222, 215, 207, 198, 189, 180, 171, 164, 158, 152, 146, 140, 134, 128, 123, 117, 111, 106, 101, 97])*10**(-3)
lnA = np.log(A)
un_sur_A = 1/A



##Affichage


plt.figure("MÃ©thode diffÃ©rentielle pour la dÃ©termination de l'ordre")
plt.clf()

ax1 = plt.axes([0.08, 0.15, 0.25, 0.7])
ax2 = plt.axes([0.415, 0.15, 0.25, 0.7])
ax3 = plt.axes([0.73, 0.15, 0.25, 0.7])

ax1.set_title("Test de l'ordre 0")
ax2.set_title("Test de l'ordre 1")
ax3.set_title("Test de l'ordre 2")

ax1.set_xlabel('$t$ (s)', fontsize = 26)
ax2.set_xlabel('$t$ (s)', fontsize = 26)
ax3.set_xlabel('$t$ (s)', fontsize = 26)

ax1.set_xticks([0, 100, 200, 300])
ax2.set_xticks([0, 100, 200, 300])
ax3.set_xticks([0, 100, 200, 300])


ax1.set_ylabel('$A$', fontsize = 26)
ax2.set_ylabel('ln($A$)', fontsize = 26)
ax3.set_ylabel('$1/A$', fontsize = 26)


l1, = ax1.plot(t, A, 'ob', markersize = 10, lw = 4)     #Trace des trois tests
l2, = ax2.plot(t, lnA, 'ob', markersize = 10, lw = 4)
l3, = ax3.plot(t, un_sur_A, 'ob', markersize = 10, lw = 4)

ax1.grid()
ax2.grid()
ax3.grid()

ax1.text(175, 0.225, '$V_{\mathrm{ClO^-}}$ = 3 mL', color = 'b', fontsize = 28)


mng = plt.get_current_fig_manager()     #Plein ecran
mng.window.showMaximized()
plt.show()


## Ajout de la regression lineaire

#fonction f realisant le fit. ici : une droite
def f(x,p) :
    (a,b) = p
    res = a*x+b
    fit = '-k*t+b'      #Legende du fit
    return (res,fit)

#derivee de la fonction f
def Dx_f(x,p) :
    (a,b) = p
    return a

#fonction d'ecart des donnees a l'ajustement
def residual(p, y, x) :
    return(y-f(x,p)[0])/np.sqrt(uy**2 + (Dx_f(x,p)*ux)**2)
    
def affichage(a, ua, k=2) :
    """Prend en arguments deux nombres et retourne deux str en ecriture scientifique a la meme puissance.
    k est le facteur de Student, vaut 2 par defaut."""
    
    dec_a = int(np.floor(np.log10(np.abs(a))))
    dec_ua = int(np.floor(np.log10(k*ua)))
    diff = int(np.abs(dec_a-dec_ua))
    
    if diff != 0 :
        if dec_a != 0 :
            str_tot = "({}$\pm${})".format(round(a*10**(-dec_a), diff), round(k*ua*10**(-dec_a),diff)) + "$\cdot 10^{"+str(dec_a)+"}$"
        else :
            str_tot = "{}$\pm${}".format(round(a*10**(-dec_a), diff), round(k*ua*10**(-dec_a),diff))
    else :
        if dec_a != 0 :
            str_tot = "({}$\pm${})".format(int(a*10**(-dec_a)), int(k*ua*10**(-dec_a))) + "$\cdot 10^{"+str(dec_a)+"}$"
        else :
            str_tot = "{}$\pm${}".format(int(a*10**(-dec_a)), int(k*ua*10**(-dec_a)))
    
    return str_tot

color_data = ['r', 'b', 'g', 'y']       #Couleurs de l'affichage
color_fit = ['c', 'orange', (0.2, 0.8, 0.5), (0.5, 0.5, 0.1)]
mark = ['o', 's', '>','<']
ref = "Mesure"
k = 2   #definition du facteur k (coeff de Student)

x_tot = t        #Variable de mesure
y_tot = lnA         #Mesurande

ux_tot = np.ones(len(x_tot)) * 20e-3/np.sqrt(3)     #J'ai gonfle les incertitudes pour que le resultat ne soit pas trop dans l'abus 
uy_tot = np.ones(len(x_tot)) * 0.02                 #Pas celle-la, elle est coherente

x = x_tot       #Les donnees x_tot et y_tot sont celles fittees
y = y_tot


#Incertitudes-types associees (et non elargies, attention)
ux = ux_tot
uy = uy_tot

popt = np.zeros((2))
pcov = np.zeros((2, 2))
upopt = np.zeros((2))
chi2r = np.zeros((1))

#estimation initale des paremetres de l'ajustement
p0 = np.array([0, 0])


#En utilisant la methpode des moindres carres sur la fonction residuelle :
result = spo.leastsq(residual, p0, args = (y,x), full_output = True)    	#args sert a completer les arguments de la fct residual, full_output de retourner la matrice de covariances
popt = result[0]	#parametres optimaux d'ajustement
pcov = result[1]	#Matrice des covariances
upopt = np.sqrt(np.abs(np.diagonal(pcov)))  #incertitude-type des paramÃ¨tres optimisees de l'ajustement
xrange = np.linspace(0,400,100)
(yrange,lab_fit) = f(xrange, popt)    	#Trace de l'ajustement

str_a = affichage(-popt[0], upopt[0])
str_b = affichage(popt[1], upopt[1])

ax2.plot(xrange, yrange, '--', color = color_fit[0], linewidth = 6, zorder = 0, label = 'y='+lab_fit+', k=' + str_a + '[L/mol/s]')
ax2.legend(framealpha = 0.5, fontsize = 14)

##Valeurs
"""Les valeurs sont issues de Gruber : on refait tout pareil pour un autre set"""

t2 = np.arange(0,381,20)+13.4
A2 = np.array([222, 203, 182, 164, 148, 130, 116, 104, 93, 83, 74, 66, 60, 53, 47, 42, 37, 33, 30, 27])*10**(-3)
lnA2 = np.log(A2)
un_sur_A2 = 1/A2


##

l4, = ax1.plot(t2, A2, 'or', markersize = 10, lw = 4)
l5, = ax2.plot(t2, lnA2, 'or', markersize = 10, lw = 4)
l6, = ax3.plot(t2, un_sur_A2, 'or', markersize = 10, lw = 4)

ax1.text(175, 0.200, '$V_{\mathrm{ClO^-}}$ = 10 mL', color = 'r', fontsize = 28)

## Ajout de la regression lineaire




x_tot = t2        #Variable de mesure
y_tot = lnA2        #Mesurande

ux_tot = np.ones(len(x_tot)) * 20e-3/np.sqrt(3)     #J'ai gonfle les incertitudes pour que le resultat ne soit pas trop dans l'abus 
uy_tot = np.ones(len(x_tot)) * 0.02                 #Pas celle-la, elle est coherente

x = x_tot
y = y_tot


#Incertitudes-types associees (et non elargies, attention)
ux = ux_tot
uy = uy_tot

popt = np.zeros((2))
pcov = np.zeros((2, 2))
upopt = np.zeros((2))
chi2r = np.zeros((1))

#estimation initale des paremetres de l'ajustement
p0 = np.array([0, 0])


#En utilisant la methpode des moindres carres sur la fonction residuelle :
result = spo.leastsq(residual, p0, args = (y,x), full_output = True)    	#args sert a completer les arguments de la fct residual, full_output de retourner la matrice de covariances
popt = result[0]	#parametres optimaux d'ajustement
pcov = result[1]	#Matrice des covariances
upopt = np.sqrt(np.abs(np.diagonal(pcov)))  #incertitude-type des paramÃ¨tres optimisees de l'ajustement
xrange = np.linspace(0,400,100)
(yrange,lab_fit) = f(xrange, popt)    	#Trace de l'ajustement

str_a = affichage(-popt[0], upopt[0])
str_b = affichage(popt[1], upopt[1])

ax2.plot(xrange, yrange, '--', color = color_fit[1], linewidth = 6, zorder = 0, label = 'y='+lab_fit+', k=' + str_a + '[L/mol/s]')
ax2.legend(framealpha = 0.5, fontsize = 14)