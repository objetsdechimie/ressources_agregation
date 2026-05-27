import matplotlib.pyplot as plt


# DONNÉES EXISTANTES


# Courbe bulle (liquide)
x_bulle = []
T_bulle = []

# Courbe rosée (vapeur)
y_rosee = []
T_rosee = []


# NOUVEAUX POINTS EXPÉRIMENTAUX


# nouveaux points liquide
x_new = []
T_x_new = []

# nouveaux points vapeur
y_new = []
T_y_new = []


# TRACÉ


plt.figure(figsize=(9,6))

# Courbe bulle
plt.plot(
    x_bulle,
    T_bulle,
    color='blue',
    marker='o',
    label='Courbe bulle (existante)'
)

# Courbe rosée
plt.plot(
    y_rosee,
    T_rosee,
    color='green',
    marker='s',
    label='Courbe rosée (existante)'
)

# Nouveaux points liquide
plt.scatter(
    x_new,
    T_x_new,
    color='red',
    marker='x',
    s=100,
    label='Nouveaux points liquide'
)

# Nouveaux points vapeur
plt.scatter(
    y_new,
    T_y_new,
    color='purple',
    marker='D',
    s=80,
    label='Nouveaux points vapeur'
)


# MISE EN FORME


plt.xlabel("Fraction molaire en propanol ($x_2$ ou $y_2$)")
plt.ylabel("Température (°C)")

plt.title("Diagramme binaire liquide-vapeur Eau / Propanol")

plt.xlim(0, 1)
plt.grid(True)

plt.legend(loc='best')

plt.tight_layout()


plt.show()