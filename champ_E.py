### Création de cartes de champ électrostatique
### créé par une ou deux particules chargées
### (Seules la direction et le sens sont affichés)
import numpy as np
import matplotlib.pylab as plt


### Constantes du problème
xmin,ymin,xmax,ymax=-2,-2,2.125,2.125 # Taille de la zone en m
h=0.25 # Pas de la zone en m
K=8.99E9 # Constante de Coulomb en USI
e=1.60E-19 # Charge élémentaire en C


### Coordonnées (x,y) des points de calcul
XX=np.arange(xmin,xmax,h)
YY=np.arange(ymin,ymax,h)
x,y=np.meshgrid(XX,YY)


### Initialisation
Ex2=None


### Particule 1
q1=-e # Charge électrique
x1=-1 # Abscisse en m
y1=0 # Ordonnée en m
### Distance particule 1-point (x,y)
r1=((x-x1)**2+(y-y1)**2)**0.5
### Coordonnée radiale champ E1
E1=K*q1/r1**2

### Coordonnées vecteur unitaire 1
ux1=(x-x1)/r1
uy1=(y-y1)/r1
### Coordonnées cartésiennes champ E1
Ex1=E1*ux1
Ey1=E1*uy1



### Particule 2 (partie à compléter)
q2=e # Charge électrique
x2=1 # Abscisse en m
y2=0 # Ordonnée en m
### Distance particule 1-point (x,y)
r2=((x-x2)**2+(y-y2)**2)**0.5
### Coordonnée radiale champ E1
E2=K*q2/r2**2

### Coordonnées vecteur unitaire 1
ux2=(x-x2)/r2
uy2=(y-y2)/r2
### Coordonnées cartésiennes champ E1
Ex2=E2*ux2
Ey2=E2*uy2


### Coordonnées du champ E
if Ex2 is None:
 Ex=Ex1
 Ey=Ey1
else:
 Ex=Ex1+Ex2
 Ey=Ey1+Ey2


### Norme des vecteurs E
E=(Ex**2+Ey**2)**0.5


### vecteur E normalisé
### (pour avoir des flèches de taille identique)
Ex=Ex/E
Ey=Ey/E


### Positions des différentes particules sources
if Ex2 is None:
 X_particules=np.array([x1])
 Y_particules=np.array([y1])
else:
 X_particules=np.array([x1,x2])
 Y_particules=np.array([y1,y2])


### Tracé des vecteurs normalisés
plt.rcParams['axes.formatter.use_locale']=True
plt.figure(1, figsize=(5, 5))
plt.quiver(x,y,Ex,Ey)
plt.plot(X_particules,Y_particules,"r.",label="Particule(s)")
plt.legend(loc='upper right')
plt.xlabel("x (en m)")
plt.ylabel("y (en m)")
plt.title("Directions et sens du champ électrostatique")
plt.show()