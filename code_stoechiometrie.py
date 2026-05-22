import numpy as np
import matplotlib.pyplot as plt

# Constante d'équilibre (choisir une valeur typique, ex : K = 1)
K = 1e5

# Fonction xNH3(ξ, a)
def x_NH3(xi, a):
    return 2 * xi / (1 + a - 2*xi)

# Équation d'équilibre à résoudre :
# K = [ 4ξ² (1 + a - 2ξ)² ] / [ (a - 3ξ)³ (1 - ξ) ]
def equilibrium_xi(a, K, P):
    # Limite supérieure pour ξ : min(a/3, 1) avec une marge de sécurité
    xi_max = min(a/3 - 1e-6, 1 - 1e-6)
    if xi_max <= 1e-6:
        return 1e-6
    
    xi_vals = np.linspace(1e-6, xi_max, 2000)
    expr = 4*xi_vals**2 * (1 + a - 2*xi_vals)**2 / ((a - 3*xi_vals)**3 * (1 - xi_vals) * P**2)
    
    # On cherche la valeur de ξ qui satisfait l'équilibre
    diff = np.abs(expr - K)
    xi_eq = xi_vals[np.argmin(diff)]
    return xi_eq

# Balayage des rapports initiaux a
a_vals = np.linspace(1, 10, 200)  # ex : de 1 à 10
x_eq_vals = [[] for _ in range(5)]  # Initialisation correcte de la liste 2D
P_values = [1e-2, 1e-1, 1, 10, 1000]

# Calcul des fractions molaires d'équilibre
for i, P in enumerate(P_values):
    for a in a_vals:
        xi_eq = equilibrium_xi(a, K, P)
        x_eq_vals[i].append(x_NH3(xi_eq, a))

# Calcul et affichage des maximas
print("Position des maximas (valeur de 'a' correspondante) :")
print("-" * 50)

maximas_positions = []
for i, P in enumerate(P_values):
    # Trouver l'index du maximum
    max_idx = np.argmax(x_eq_vals[i])
    max_position = a_vals[max_idx]
    max_value = x_eq_vals[i][max_idx]
    maximas_positions.append(max_position)
    
    print(f"P = {P:>8} : a_max = {max_position:.3f}, x_NH3_max = {max_value:.6f}")

print("-" * 50)

# Tracé des courbes
plt.figure(figsize=(10, 6))
colors = ['blue', 'red', 'green', 'orange', 'purple']
labels = [f'P = {P}' for P in P_values]

for i in range(len(P_values)):
    plt.plot(a_vals, x_eq_vals[i], color=colors[i], label=labels[i], linewidth=2)
    
    # Marquer le maximum sur la courbe
    max_idx = np.argmax(x_eq_vals[i])
    plt.plot(a_vals[max_idx], x_eq_vals[i][max_idx], 
             marker='o', color=colors[i], markersize=8, markerfacecolor='white', 
             markeredgecolor=colors[i], markeredgewidth=2)
    
    # Annoter le maximum
    plt.annotate(f'a={a_vals[max_idx]:.2f}', 
                xy=(a_vals[max_idx], x_eq_vals[i][max_idx]),
                xytext=(5, 5), textcoords='offset points',
                fontsize=14, color=colors[i], fontweight='bold')

plt.xlabel(r"Rapport initial $a = \dfrac{n_{H_2}^0}{n_{N_2}^0}$")
plt.ylabel(r"Fraction molaire $x_{NH_3}$ à l'équilibre")
plt.title("Évolution de $x_{NH_3}$ en fonction du rapport initial des réactifs")
plt.legend()
plt.grid(True, alpha=0.3)
plt.tight_layout()
plt.show()