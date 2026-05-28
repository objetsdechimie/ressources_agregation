#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#%matplotlib qt  
"""
Python code provided as is.
Made by Thibault Fogeron, and based on code developped by Vincent Wieczny from Chemistry Department, ENS de Lyon, France
This code is under licence CC-BY-NC-SA. It enables you to reuse the code by mentioning the orginal author and without making profit from it.
"""

#Librairies
import matplotlib.pyplot as plt
import numpy as np
import widgets
import scipy.constants as constants
from matplotlib import rc

################################
### Paramater initialization ###
################################

#Physical constants
F=96500.0 #Faraday number (C/mol)
R=8.314 #Gas constant (J/K/mol)
T=298.0 #Temperature (K)


#Redox system data

#Ce4+/Ce3+
E_std_Ce=1.44 #standard potential (V/SHE)
D_Ce3=1e-9
D_Ce4=1e-9
c0_Ce4=1e-3 #titrant concentration (mol/L)

#Fe3+/Fe2+
E_std_Fe=0.77 #standard potential (V/SHE)
D_Fe3=1.5e-9
D_Fe2=1.5e-9
c0_Fe2=1e-3 #titrated concentration (mol/L)

#Electrochemical setup and system
delta=1e-5 #diffuse layer thickness (m)
n=1 #number of exchanged electrons
A=1e-5 #electrode area (m²)
alpha = 0.5 #Buttler-Volmer coefficient

#Titration parameters
V_Fe2=10 #titrated volume (mL)
V0=10 #total volume (mL)

## Modulated parameters

parameters = {'V' : widgets.FloatSlider(value=0.000001, description='$V$ $\mathrm{(mL)}$', min=0.00000001, max=20),
              'k1' : widgets.IntSlider(value=1, description='$k_{Fe}$', min=0, max=1),
              'k2' : widgets.IntSlider(value=1, description='$k_{Ce}$', min=0, max=1),
              'E0' : widgets.FloatSlider(value=1.000, description='$E_0$ $\mathrm{(V/ESH)}$', min=0.00000001, max=2)
              }

# V is the added volume during titration
# k1 and k2 are on/off switch that allow to describe either fast or slow charge kinteic transferts
# E0 is the value of the ptentiel imopsed

#################
### Functions ###
#################

# Titration volume determination
def Veq(c0_Fe2,c0_Ce4,V_Fe2):
    return c0_Fe2*V_Fe2/c0_Ce4

# Fe2+ concentration as a function of the added volume V
def c_Fe2(V,Veq):
    if V<=Veq:
        return (c0_Fe2*V_Fe2-c0_Ce4*V)/(V0)
    else:
        return 0.0000000001 # non-zero value due to a divergent behaviour

# Fe3+ concentration as a function of the added volume V
def c_Fe3(V,Veq):
    if V<Veq:
        return c0_Ce4*V/(V0)
    else:
        return c0_Ce4*Veq/(V0)

# Ce4+ concentration as a function of the added volume V
def c_Ce4(V,Veq):
    if V<Veq:
        return 0.0000000001
    else:
        return c0_Ce4*(V-Veq)/(V0)

# Ce3+ concentration as a function of the added volume V
def c_Ce3(V,Veq):
    if V<Veq:
        return c0_Ce4*V/(V0)
    else:
        return c0_Ce4*Veq/(V0)

#Anodic diffusive controlled current
def i_a(delta,Cred,Dred):
  return n*F*A*Dred*Cred/delta

#Cathodic diffusive controlled current
def i_c(delta,Cox,Dox):
  return -n*F*A*Dox*Cox/delta

#Diffusion-limited current
def i_diff(delta,Dred,Dox,Cred,Cox,E,E_std):
    ia=i_a(delta,Cred,Dred)
    ic=i_c(delta,Cox,Dox)
    k=(E-E_std)*(n*F)/(R*T)
    return (np.exp(k)*ia+ic)/(1+np.exp(k))

#Electronic transfer-limited current


def kanodicreduced(E_std,E,alpha,F,R,T):
    """
    anodic kinetic constant without the k0 prefactor see 7.16 of Girault
    """
    NernstCoeff =  F / (R * T)
    return np.exp(alpha * (E-E_std)*NernstCoeff)

def kcathodicreduced(E_std,E,alpha,F,R,T):
    """
    cathodic kinetic constant without the k0 prefactor see 7.16 of Girault
    """
    NernstCoeff =  F / (R * T)
    return np.exp(-(1-alpha) * (E-E_std)*NernstCoeff)



def currentButlerVolmer(E_std, alpha, E,F,R,T,Cred,Cox,A):


    fk0 = 10**(-7) #Value of k0 arbitrary
    ka = kanodicreduced(E_std,E,alpha,F,R,T)
    kc = kcathodicreduced(E_std,E,alpha,F,R,T)
    i = n * F * A * fk0 * (Cred * ka -Cox* kc )
    return i

#i-E Fe data at the added volume V:
def iE_data_Fe(V,Veq,E,k1):
    #Concentration calculation
    Fe2=c_Fe2(V,Veq)
    Fe3=c_Fe3(V,Veq)

# k1 is just a switch to go from fast to slow electron transfert
#fast kinetic for the electron transfert : Diffusion-limited current only
    if k1 >=1:
        i_diff_Fe=i_diff(delta,D_Fe2,D_Fe3,Fe2,Fe3,E,E_std_Fe)
        return i_diff_Fe

#slow kinetic for the electron transfert : average of Diffusion-limited current and Buttler Volmer

    if k1 < 1:
        i_diff_Fe=i_diff(delta,D_Fe2,D_Fe3,Fe2,Fe3,E,E_std_Fe)
        i_butt_Fe=currentButlerVolmer(E_std_Fe, alpha, E,F,R,T,Fe2,Fe3,A)
        i_Fe = i_diff_Fe*i_butt_Fe/(i_diff_Fe+i_butt_Fe)
        return i_Fe


#i-E data at the added volume V:
def iE_data_Ce(V,Veq,E,k2):
    #Concentration calculation
    Ce4=c_Ce4(V,Veq)
    Ce3=c_Ce3(V,Veq)

#fast kinetic for the electron transfert : Diffusion-limited current only
    if k2 >=1:
        i_diff_Ce=i_diff(delta,D_Ce3,D_Ce4,Ce3,Ce4,E,E_std_Ce)
        return i_diff_Ce

#slow kinetic for the electron transfert : average of Diffusion-limited current and Buttler Volmer

    if k2 < 1:
        i_diff_Ce=i_diff(delta,D_Ce3,D_Ce4,Ce3,Ce4,E,E_std_Ce)
        i_butt_Ce=currentButlerVolmer(E_std_Ce, alpha, E,F,R,T,Ce3,Ce4,A)
        i_Ce = i_diff_Ce*i_butt_Ce/(i_diff_Ce+i_butt_Ce)
        return i_Ce



#i-E data at the added volume V:
def iE_data_tot(V,Veq,E,k1,k2):

    #Curve calculation
    i_Fe=iE_data_Fe(V,Veq,E,k1)
    i_Ce=iE_data_Ce(V,Veq,E,k2)
    i_tot=i_Fe+i_Ce

    return i_tot



#Titration spot


def titration_i(V,E0,k1,k2):
    counter=0
    for k in range(0,len(E)):
        if E[k]<E0:
            counter=counter+1
    return iE_data_tot(V,Veq,E,k1,k2)[counter]

# E which is defined as a list is the absciss of the i= f(E) curves
# the function for allow to convert the potential E0 in one of the element of the list
# then we calculate the corresponding value of i=f(E0)

#Titration curve
def titration_curve(V_domain,E0,k1,k2):
    i_titr=[]
    for v in V_domain:
        i_titr.append(titration_i(v,E0,k1,k2))
    return i_titr


#===========================================================
# --- Initialization of the plot ---------------------------
#===========================================================

#fig,(ax1,ax2)=plt.subplots(1,2,figsize=(16,6))
fig=plt.figure(figsize=(18,6))

fig.suptitle(r'Amperometric titration $(E=E_0)$  of a $\mathbf{Fe^{2+}}$ solution by a $\mathbf{Ce^{4+}}$ solution',weight='bold')


Veq=Veq(c0_Fe2,c0_Ce4,V_Fe2)


E=np.arange(0.0001,2.001,0.0211)
V_domain=np.arange(0.0000001,2*Veq+0.00001,0.02)



fig.text(0.01,0.9,r'Titration conditions', multialignment='left', verticalalignment='top',weight='bold')
fig.text(0.01,0.85,r'Titrated solution', multialignment='left', verticalalignment='top')
fig.text(0.01,0.80,r'$c_0=${:.3f} mol/L'.format(c0_Fe2), multialignment='left', verticalalignment='top')
fig.text(0.01,0.77,r'$V_0=${:.2f} mL'.format(V_Fe2), multialignment='left', verticalalignment='top')
fig.text(0.01,0.72,r'Titrant solution', multialignment='left', verticalalignment='top')
fig.text(0.01,0.67,r'$c_1=${:.3f} mol/L'.format(c0_Ce4), multialignment='left', verticalalignment='top')
fig.text(0.01,0.62,r'End point', multialignment='left', verticalalignment='top')
fig.text(0.01,0.57,r'$V_e=${:.2f} mL'.format(Veq), multialignment='left', verticalalignment='top')
#fig.text(0.01,0.52,r'$i_0=${:.3f} µA'.format(i0bis), multialignment='left', verticalalignment='top')

fig.text(0.01,0.47,r'Dilution is not taken into account.', multialignment='left', verticalalignment='top')
fig.text(0.60,0.09,r'$k_i = 1 $ => fast electronic transfert', multialignment='left', verticalalignment='top')
fig.text(0.60,0.07,r'$k_i = 0 $ => slow electronic transfert', multialignment='left', verticalalignment='top')


ax1 = fig.add_axes([0.2, 0.2, 0.35, 0.7])
ax2 = fig.add_axes([0.60, 0.2, 0.35, 0.7])
#ax.axhline(0, color='k')



ax1.text(0.77,2.25e-7,'$\mathrm{Fe^{3+}_{(aq)} + e^- \leftrightarrows \, Fe^{2+}_{(aq)}}$',horizontalalignment='center',
     verticalalignment='center')

ax1.text(1.44,2.25e-7,'$\mathrm{Ce^{4+}_{(aq)} + e^- \leftrightarrows \, Ce^{3+}_{(aq)}}$',horizontalalignment='center',
     verticalalignment='center')

ax1.plot([E.min(), E.max()],[0,0],':',lw=1,color='grey')

ax1.set_xlim(E.min(), E.max())
ax1.set_ylim(-3e-7,3e-7)

ax1.set_xlabel('$E$ $\mathrm{(V/ESH)}$')
ax1.set_ylabel('Current $i$ $\mathrm{(A)}$')

ax2.plot([Veq,Veq],[-3e-7,3e-7],':',lw=1,color='grey')

ax2.plot([V_domain.min(), V_domain.max()],[0,0],':',lw=1,color='grey')



ax2.set_xlim(V_domain.min(), V_domain.max())
ax2.set_ylim(-3e-7,3e-7)


ax2.set_xlabel('$V$ $\mathrm{(mL)}$')
ax2.set_ylabel('$i$ $\mathrm{(A)}$')


#===========================================================
# --- Plot of the updated curves ---------------------------
#===========================================================

# This function is called when the sliders are changed
def plot_data(V,E0,k1,k2):
    lines['$i_\mathrm{Fe}$'].set_data(E,iE_data_Fe(V,Veq,E,k1))
    lines['$i_\mathrm{Ce}$'].set_data(E,iE_data_Ce(V,Veq,E,k2))
    lines['$i_\mathrm{tot}$'].set_data(E,iE_data_tot(V,Veq,E,k1,k2))
    lines['$E_0$'].set_data([E0,E0],[-3e-7,3e-7])
    lines['$Titration \ curve$'].set_data(V_domain,titration_curve(V_domain,E0,k1,k2))
    lines['$Titration \ spot \ (left)$'].set_data(E0,titration_i(V,E0,k1,k2))
    lines['$Titration \ spot \ (right)$'].set_data(V,titration_i(V,E0,k1,k2))

    fig.canvas.draw_idle()

# See below for the description of the plot






lines = {}
#Plot of the i =f(E) curves (left)
lines['$i_\mathrm{Fe}$'], = ax1.plot([], [],color='green',lw=2,label='$i_\mathrm{Fe}$')
lines['$i_\mathrm{Ce}$'], = ax1.plot([], [],color='blue',lw=2,label='$i_\mathrm{Ce}$')
lines['$i_\mathrm{tot}$'], = ax1.plot([], [], lw=3, color='red',label='$i_\mathrm{tot}$')
#Display of the i(E0) value
lines['$E_0$'],=ax1.plot([], [],'--',color='grey',lw=1,label=r'$E_0$')
lines['$Titration \ spot \ (left)$'], = ax1.plot([], [],'o',color='black',lw=2,label='$Titration \ spot$')
#Titration curve and evolution of the i(E0) value during titration
lines['$Titration \ curve$'], = ax2.plot([], [],color='red',lw=3,label='$Titration \ curve$')
lines['$Titration \ spot \ (right)$'], = ax2.plot([], [],'o',color='black',lw=2,label='$Titration \ spot$')


ax1.legend()


ax2.legend()



param_widgets = widgets.make_param_widgets(parameters, plot_data, slider_box=[0.20, 0.02, 0.35, 0.1])
choose_widget = widgets.make_choose_plot(lines, box=[0.01,0.2,0.12, 0.2])
reset_button = widgets.make_reset_button(param_widgets,box=[0.85, 0.05, 0.10, 0.05])

if __name__=='__main__':
    plt.show()