import numpy as np


# ETALONNAGE - 3 PRODUITS


produits = ["Produit 1", "Produit 2", "Produit 3"]

# Aires analytes étalon
A_std = np.array([
    154230,
    132500,
    175800
])

# Aires étalon interne étalon
A_ei_std = np.array([
    132540,
    132540,
    132540
])

# Concentrations analytes étalon (mol/L)
C_std = np.array([
    25.0,
    18.0,
    30.0
])

# Concentration étalon interne (mol/L)
C_ei_std = np.array([
    20.0,
    20.0,
    20.0
])


# INCERTITUDES ETALON


u_A_std = np.array([
    1200,
    1100,
    1300
])

u_A_ei_std = np.array([
    1000,
    1000,
    1000
])

u_C_std = np.array([
    0.20,
    0.15,
    0.25
])

u_C_ei_std = np.array([
    0.15,
    0.15,
    0.15
])


# CALCUL FACTEURS DE REPONSE


F = (A_std * C_ei_std) / (A_ei_std * C_std)

u_rel_F = np.sqrt(
    (u_A_std / A_std)**2 +
    (u_A_ei_std / A_ei_std)**2 +
    (u_C_std / C_std)**2 +
    (u_C_ei_std / C_ei_std)**2
)

u_F = F * u_rel_F


# DONNEES ECHANTILLONS


# Aires analytes échantillon
A_x = np.array([
    184520,
    145800,
    210450
])

# Aires étalon interne échantillon
A_ei = np.array([
    143210,
    143210,
    143210
])

# Concentration EI échantillon (mol/L)
C_ei = np.array([
    20.0,
    20.0,
    20.0
])


# INCERTITUDES ECHANTILLON


u_A_x = np.array([
    1500,
    1200,
    1600
])

u_A_ei = np.array([
    1100,
    1100,
    1100
])

u_C_ei = np.array([
    0.15,
    0.15,
    0.15
])


# CALCUL CONCENTRATIONS


C_x = (A_x / A_ei) * (C_ei / F)


# PROPAGATION INCERTITUDES


u_rel_Cx = np.sqrt(
    (u_A_x / A_x)**2 +
    (u_A_ei / A_ei)**2 +
    (u_C_ei / C_ei)**2 +
    (u_F / F)**2
)

u_Cx = C_x * u_rel_Cx

# Incertitude élargie
U_Cx = 2 * u_Cx


# AFFICHAGE

#les :.3f permettent de régler le nombre de chiffres significatifs !!

for i in range(3):

    print(f"\n{produits[i]}")

    print(f"Facteur de réponse F = {F[i]:.5f}")
    print(f"u(F) = {u_F[i]:.5f}")

    print(f"Concentration = {C_x[i]:.3f} mol/L")

    print(f"u(Cx) = {u_Cx[i]:.3f} mol/L")

    print(f"U(Cx) élargie (k=2) = ±{U_Cx[i]:.3f} mol/L")

