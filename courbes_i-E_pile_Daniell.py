import numpy as np
import matplotlib.pyplot as plt

"""
Tracé des courbes intensité-potentiel pour la pile Daniell
"""

#%% Définition des constantes

i0 = 0.01
i0_eau = 0.00001
alpha = 0.5
n = 2 #nombre d'électrons échangés
F = 96500
R = 8.314
T = 298
A = 1e-4
DO=6e-6
DR=6e-6
delta=1e-6
cR=0.01
cO=0.01
E_Cu=0.34 #potentiel standard du couple Cu2+/Cu [en V]
E_Zn=-0.76 #potentiel standard du couple Zn2+/Zn [en V]
eta_eau=-0.90 #surtension du couple H+/H2 sur une électrode de zinc [en V] On considère que c'est la même sur le cuivre.

f = F/(R*T)
E_cat=E_Cu+R*T/(n*F)*np.log(cO)
E_an=E_Zn+R*T/(n*F)*np.log(cR)
mO=DO/delta
mR=DR/delta

#%% Définition de l'abscisse et des ordonnées

E0 = np.linspace(-1.5,eta_eau,100)
E1 = np.linspace(eta_eau,E_an,10)
E2 = np.linspace(E_an,E_cat,100)
E3 = np.linspace(E_cat,0.75,100)
E=np.concatenate((E0,E1,E2,E3))

N=len(E)
abscisses=np.zeros(N)


#courants limites de diffusion
ial = F*A*mR*cR
icl = -F*A*mO*cO

#Pour les courants anodique/cathodique/totaux en transfert de charge limitant

def i_charge(i0,E,Es):
    #Contrôle de charges limitant
    return i0*np.exp(alpha*n*f*(E-Es))-i0*np.exp(-(1-alpha)*n*f*(E-Es))

def i_diff(i0,E,Es):
    #Pour la diffusion limitante
    return (ial*np.exp(n*f*(E-Es)) + icl)/(1+np.exp(n*f*(E-Es)))

def it(i0,E,Es):
    #Pour un contrôle mixte
    np.seterr(divide='ignore', invalid='ignore')
    return (i_charge(i0,E,Es)*i_diff(i0,E,Es))/(i_charge(i0,E,Es)+i_diff(i0,E,Es))

#Pour l'anode :
i_an0=i_charge(i0_eau,E0,eta_eau)+it(i0,E0,E_an)
i_an1=it(i0,E1,E_an)
i_an2=i_charge(i0,E2,E_an)

#Pour la cathode :
i_cat0=i_charge(i0_eau,E0,eta_eau)+it(i0,eta_eau,E_cat)
i_cat1=it(i0,E1,E_cat)
i_cat2=it(i0,E2,E_cat)
i_cat3=i_charge(i0,E3,E_cat)


#%% Tracé du graphique

plt.grid()
plt.title("Courbes intensité-potentiel pour la pile Daniell")
plt.xlabel("Potentiel E (V)")
plt.ylabel("Intensité du courant i (A)")
plt.xlim(-1.5,0.75)
plt.ylim(-0.8,1)
plt.text(E_an,-0.1,'$E_{an}$ = -0,82 V', horizontalalignment = 'center', verticalalignment = 'center')
plt.text(E_cat,-0.1,'$E_{cat}$ = 0,28 V', horizontalalignment = 'center', verticalalignment = 'center')
plt.text(-0.68,0.5,'Zn -> Zn$^{2+}$', horizontalalignment = 'center', verticalalignment = 'center')
plt.text(0.2,-0.3,'Cu $\leftarrow$ Cu$^{2+}$', horizontalalignment = 'center', verticalalignment = 'center')
plt.text(-1.13,-0.7,'H$_2\leftarrow$ H$^{+}$', horizontalalignment = 'center', verticalalignment = 'center')

#Axe des abscisses
plt.plot(E,abscisses,color='k')

#Pour l'anode :
plt.plot(E0,i_an0,'g--')
plt.plot(E1,i_an1,'g--')
plt.plot(E2,i_an2,color="g",label="à l'anode de zinc")

#Pour la cathode :
plt.plot(E0,i_cat0,'b',label="à la cathode de cuivre")
plt.plot(E1,i_cat1,'b')
plt.plot(E2,i_cat2,'b')
plt.plot(E3,i_cat3,"b--")

plt.legend()
plt.show()
