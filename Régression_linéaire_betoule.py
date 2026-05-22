import numpy as np
import matplotlib.pyplot as plt
 
# =============================================================================
#  INTERFACE UTILISATEUR — À MODIFIER SELON VOTRE SITUATION EXPÉRIMENTALE
# =============================================================================
 
TITRE = "Détermination de G"      # Titre de l'expérience
X_LABEL = "w^2 r^2"  # Légende axe X (avec unité)
Y_LABEL = "z"              # Légende axe Y (avec unité)
 
# --- Séries de données (même longueur obligatoire) ---
X  = [0.06516, 0.1080, 0.1568, 0.2073, 0.2710]
sX = [0.00296,0.00382,0.00460,0.00529,0.00605]  # Incertitudes-types sur X
 
Y  = [0.0040,0.0068,0.0093,0.0127,0.0156]
sY = [0.001, 0.001, 0.001, 0.001, 0.001]                         # Incertitudes-types sur Y
 
# --- Paramètres de la simulation ---
N_SIMULATIONS = 10000   # Nombre de tirages Monte-Carlo (recommandé : 10 000)
T_STUDENT     = 2.571   # Coefficient de Student pour IC à 95%
                        # Exemples courants (alpha=0.05) :
                        #   ν=3  (5 pts)  -> 3.182
                        #   ν=4  (6 pts)  -> 2.776
                        #   ν=5  (7 pts)  -> 2.571
                        #   ν=6  (8 pts)  -> 2.447
                        #   ν=8  (10 pts) -> 2.306
                        #   ν=inf         -> 1.960
 
# =============================================================================
#  CALCULS — NE PAS MODIFIER
# =============================================================================
 
X  = np.array(X,  dtype=float)
sX = np.array(sX, dtype=float)
Y  = np.array(Y,  dtype=float)
sY = np.array(sY, dtype=float)
 
# --- Régression linéaire par moindres carrés ordinaires ---
a, b = np.polyfit(X, Y, 1)
 
# --- Pondération khi-deux : p_i = 1 / sqrt(sY_i² + (a·sX_i)²) ---
p = 1.0 / np.sqrt(sY**2 + (a * sX)**2)
 
# --- Régression linéaire pondérée (méthode du khi-deux) ---
a_khi2, b_khi2 = np.polyfit(X, Y, 1, w=p)
 
# --- Méthode de Monte-Carlo ---
rng   = np.random.default_rng()   # Générateur reproductible : rng = np.random.default_rng(42)
Y_rd  = rng.normal(Y,  sY,  (N_SIMULATIONS, len(Y)))
X_rd  = rng.normal(X,  sX,  (N_SIMULATIONS, len(X)))
 
a_rd, b_rd = [], []
for i in range(N_SIMULATIONS):
    # Recalcul de la pondération sur chaque jeu simulé
    p_i = 1.0 / np.sqrt(sY**2 + (a_khi2 * sX)**2)
    coeffs = np.polyfit(X_rd[i], Y_rd[i], 1, w=p_i)
    a_rd.append(coeffs[0])
    b_rd.append(coeffs[1])
 
a_rd = np.array(a_rd)
b_rd = np.array(b_rd)
 
# --- Incertitudes-types et intervalles de confiance à 95 % ---
s_a = np.std(a_rd)
s_b = np.std(b_rd)
U_a = T_STUDENT * s_a
U_b = T_STUDENT * s_b
 
# --- Coefficient de détermination R² (ajustement khi-deux) ---
y_pred = a_khi2 * X + b_khi2
ss_res = np.sum((Y - y_pred)**2)
ss_tot = np.sum((Y - np.mean(Y))**2)
R2 = 1 - ss_res / ss_tot
 
# =============================================================================
#  AFFICHAGE DES RÉSULTATS
# =============================================================================
 
print("=" * 55)
print(f"  {TITRE}")
print("=" * 55)
print(f"  Équation : Y = a·X + b")
print(f"  a      = {a_khi2:.4E}")
print(f"  U(a)   = {U_a:.2E}   (IC 95 %, t_s = {T_STUDENT})")
print(f"  b      = {b_khi2:.4E}")
print(f"  U(b)   = {U_b:.2E}   (IC 95 %)")
print(f"  R²     = {R2:.6f}")
print(f"  N_sim  = {N_SIMULATIONS}")
print("=" * 55)
 
# =============================================================================
#  FIGURES
# =============================================================================
 
fig, axes = plt.subplots(1, 3, figsize=(16, 5))
fig.suptitle(TITRE, fontsize=14, fontweight="bold")
 
# --- Figure 1 : Régression avec barres d'incertitude ---
ax = axes[0]
ax.errorbar(X, Y, xerr=sX, yerr=sY,
            fmt='o', color='steelblue', ecolor='grey',
            capsize=4, linewidth=1.2, label='Données ± u')
x_fit = np.linspace(X.min() - 0.05*(X.max()-X.min()),
                    X.max() + 0.05*(X.max()-X.min()), 300)
ax.plot(x_fit, a_khi2 * x_fit + b_khi2,
        color='tomato', linewidth=2,
        label=f'Y = {a_khi2:.3E}·X + {b_khi2:.3E}\nR² = {R2:.4f}')
ax.set_xlabel(X_LABEL, fontsize=12)
ax.set_ylabel(Y_LABEL, fontsize=12)
ax.set_title("Régression pondérée (χ²)", fontsize=12)
ax.legend(fontsize=10)
ax.tick_params(labelsize=11)
ax.grid(True, linestyle='--', alpha=0.4)
 
# --- Figure 2 : Histogramme de a ---
ax = axes[1]
ax.hist(a_rd, bins=80, color='steelblue', alpha=0.75, density=True)
ax.axvline(a_khi2,         color='tomato',    linewidth=2,   label=f'a = {a_khi2:.3E}')
ax.axvline(a_khi2 - U_a,  color='darkorange', linewidth=1.5, linestyle='--', label=f'± U(a) = {U_a:.2E}')
ax.axvline(a_khi2 + U_a,  color='darkorange', linewidth=1.5, linestyle='--')
ax.set_xlabel("Coefficient a", fontsize=12)
ax.set_ylabel("Densité", fontsize=12)
ax.set_title("Distribution Monte-Carlo de a", fontsize=12)
ax.legend(fontsize=10)
ax.tick_params(labelsize=11)
ax.grid(True, linestyle='--', alpha=0.4)
 
# --- Figure 3 : Histogramme de b ---
ax = axes[2]
ax.hist(b_rd, bins=80, color='seagreen', alpha=0.75, density=True)
ax.axvline(b_khi2,         color='tomato',    linewidth=2,   label=f'b = {b_khi2:.3E}')
ax.axvline(b_khi2 - U_b,  color='darkorange', linewidth=1.5, linestyle='--', label=f'± U(b) = {U_b:.2E}')
ax.axvline(b_khi2 + U_b,  color='darkorange', linewidth=1.5, linestyle='--')
ax.set_xlabel("Ordonnée à l'origine b", fontsize=12)
ax.set_ylabel("Densité", fontsize=12)
ax.set_title("Distribution Monte-Carlo de b", fontsize=12)
ax.legend(fontsize=10)
ax.tick_params(labelsize=11)
ax.grid(True, linestyle='--', alpha=0.4)
 
plt.tight_layout()
plt.show()
 