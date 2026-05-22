# -*- coding: utf-8 -*-
"""
Created on Thu Mar 26 09:51:21 2026

@author: maxim
"""

#Michelson en lame d'air
#!!!Commentaires indiquant une variable pouvant Ãªtre modifiee!!!

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button # outils pour faire des curseurs


# setup some generic data
N = 80 # Nombre de pixels
x, y = np.mgrid[-N:N, -N:N]*10**(-2)#x=array : [[-N...-N]...[N...N]] et y = array : [[-N...N]...[-N...N]] (en mm)
I0=1#!!!Intensite de la source!!!
Lambda = 589*10**(-9) #!!!longueur d'onde en m!!!
e_init=0 #!!!decalage (en microm) entre les deux miroirs!!!
f = 0.5 #!!!Distance focale de la lentille utilisee pour projeter les interferences!!!
r2 = x**2+y**2
I = I0/2*(1+np.cos(4*np.pi/Lambda*e_init*10**(-6)*(1-r2/(2*f**2))))


fig, ax = plt.subplots()
image = ax.imshow(I, cmap='copper',vmin=0,vmax=1, origin='lower', interpolation='none')#Transforme la matrice I en une image
fig.colorbar(image, ax=ax, label='Intensite (arb. units)')#legende (axe en couleur)
ax.set_xlabel('x (en cm)')
ax.set_ylabel("y (en cm)")
ax.set_title("Michelson en lame d'air")
grad = [0,30, N, 130, 2*N]
plt.yticks(grad, ['','-50', '0' ,'50',''])
plt.xticks(grad, ['','-50', 0 , '50', ''])


# Ajuste le graphique pour faire de la place pour le curseur
plt.subplots_adjust(left=0.25, bottom=0.25)


# Creation du curseur pour controler e
ax_e = plt.axes([0.25, 0.1, 0.65, 0.03])
e_slider = Slider(  
    ax=ax_e,
    label='e ($\mu$m)',
    valmin=0,
    valmax=1,
    valinit=e_init,
)


# Fonction appelee des que la valeur de x change
def update(val):
    I = I0/2*(1+np.cos(4*np.pi/Lambda*e_slider.val*10**(-6)*(1-r2/(2*f**2))))
    image.set_data = ax.imshow(I, cmap='copper',vmin=0,vmax=1, origin='lower', interpolation='none')
    fig.canvas.draw_idle()


e_slider.on_changed(update)

# Creation d'un bouton reset
resetax = plt.axes([0.8, 0.025, 0.1, 0.04])
button = Button(resetax, 'Reset', hovercolor='0.975')


def reset(event):
    e_slider.reset()
button.on_clicked(reset)



plt.show()