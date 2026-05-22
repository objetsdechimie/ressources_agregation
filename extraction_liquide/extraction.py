# -*- coding: utf-8 -*-
"""
Created on Wed May 20 13:16:09 2026

@author: maxim
"""


"""
Ce code a pour objectif de simuler des extractions liquide-liquide et de voir les effets des diffÃ©rentes paramÃ¨tres que sont :
    - n le nombre d'extractions
    - K la constante de partage
    - V le volume d'extraction utilisÃ© 
Ce code a Ã©tÃ© Ã©crit par Raphael Rullan, citer si utilisation. Je remercie Vincent Wieczny <3 pour le module Widgets qui permet de faire les sliders. 
On s'intÃ©resse Ã  la constante de partage qui est dÃ©finie comme : K = [BOH]_{Et2O}/[BOH]_{H20}. On extraie avec l'ether diethylique. On note q la fraction restÃ©e dans l'eau'

"""

#Importation des librairies

import numpy as np
import matplotlib.pyplot as plt
import widgets

#DÃ©finition des paramÃ¨tres des sliders

parameters = {'K' : widgets.IntSlider(value=3, description='Constante de partage', min=0, max=100),
              'V_1'  : widgets.IntSlider(value = 10, description="Volume de solvant d'extraction (mL)", min = 1, max =100),
              'V_2'  : widgets.IntSlider(value = 10, description="Volume de solvant Ã  extraire (mL)", min = 1, max =100),
              'N'  : widgets.IntSlider(value = 2, description="nombre d'extraction", min = 1, max =100)}
              

#DÃ©finition des fonctions

def coefficient_extraction_n(K,V_1,V_2,n):
    '''ce qui nous intÃ©resse ici ce n'est pas la partie qui reste dans l'eau mais celle qui est extraite'''
    q  = (V_2/(V_2 + K*V_1))**n
    return 1 - q

def coefficient_extraction_V(K,x):
    q = 1/(1 + K*x)
    return 1 - q

def coefficient_extraction_unique(K,V_1,V_2,n):
    q = (V_2/(V_2 + K*V_1))**(1+n)
    return 1 - q
    
#RÃ©alisation 


n = np.linspace(1, 100, 100)
X = np.linspace(0,10,1000)
n0 = np.zeros(100)

fig=plt.figure(figsize=(12,10))
plt.subplots_adjust(left=0.125,bottom=0.2,right=0.8,top=0.9,wspace=0.3,hspace=0.5)
ax1 = plt.subplot(211)
ax2 = plt.subplot(212)

ax2.set_ylim(0,1.2)
ax1.set_xlim(0,10)
ax1.set_ylim(0.1)

def plot_data(K,V_1,V_2,N):
    lines['coeff_extract_n'].set_data(n,coefficient_extraction_n(K, V_1, V_2, n))
    lines['coeff_extract_V'].set_data(X,coefficient_extraction_V(K, X))
    lines['coeff_unique'].set_data(n,coefficient_extraction_unique(K,V_1,V_2,n0) )
    ax2.set_xlim(1,N)
    texts.set_text('r = {:.2e}'.format(coefficient_extraction_unique(K,V_1,V_2,n0)[0]))
    
lines={}
lines['coeff_extract_n'], = ax2.plot([], [],'r-',lw=2,label="extraction en fonction du nombre d'extractions")
lines['coeff_unique'], = ax2.plot([], [],'go',label="coefficient d'extraction Ã  chaque Ã©tape")
lines['coeff_extract_V'], = ax1.plot([], [],'b-',lw=2,label='extraction en fonction du rapport des volumes')

ax1.set_xlabel("rapport des volumes x")
ax2.set_xlabel("nombre d'extractions")
texts=  ax2.text(0.01, 0.05,"coefficient d'extraction =" ,bbox=dict(facecolor='red', alpha=0.5), fontsize=10, transform=ax2.transAxes)
ax2.set_ylabel("pourcentage d'espÃ¨ce extraite")
ax1.set_ylabel("pourcentage d'espÃ¨ce extraite")

ax1.set_title("Influence du rapport des volumes sur le coefficient d'extraction")
ax2.set_title("Influence du nombre d'extractions sur le coefficient d'extraction")
ax1.legend()
ax2.legend(bbox_to_anchor=(0.05,0.1, 0.8,1), loc = 4)
param_widgets = widgets.make_param_widgets(parameters, plot_data, slider_box=[0.35, 0.02, 0.35, 0.1])
choose_widget = widgets.make_choose_plot(lines, box=[0.85,0.2,0.15, 0.2])
reset_button = widgets.make_reset_button(param_widgets,box=[0.01, 0.01, 0.08, 0.05])













if __name__=='__main__':
    plt.show()