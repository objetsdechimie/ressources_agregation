
""" 

script python permettant de faire des regressions linéaires inspiré du code issu du livre Outils Numériquesa
au service de la chimie expérimentale de S.Betoule
"""



import matplotlib.pyplot as plt
import numpy as np

# =============================================================================
# CONTEXTE EXPÉRIMENTAL (PARAMÈTRES MODIFIABLES)
# =============================================================================

# Définition des étiquettes et du titre pour les rapports et graphiques
TITRE_MANIP = "Détermination de G"
NOM_AXE_X = "w^2 r^2"
NOM_AXE_Y = "z"

# Séries de mesures issues de l'expérience (les 4 listes doivent avoir la même taille)
donnees_x = [0.06516, 0.1080, 0.1568, 0.2073, 0.2710]
incertitudes_x = [
    0.00296,
    0.00382,
    0.00460,
    0.00529,
    0.00605,
]  # Écarts-types u(x)
donnees_y = [0.0040, 0.0068, 0.0093, 0.0127, 0.0156]
incertitudes_y = [0.001, 0.001, 0.001, 0.001, 0.001]  # Écarts-types u(y)

# Paramétrage de l'évaluation statistique
TOTAL_TIRAGES = 10000  # Nombre de réplications pour la méthode Monte-Carlo
FACTEUR_STUDENT = (
    2.571  # Coefficient correctif pour un intervalle de confiance à 95%
)

# =============================================================================
# PHASE DE CALCULS ET TRAITEMENT NUMÉRIQUE
# =============================================================================

# Conversion des listes en vecteurs NumPy pour activer les opérations mathématiques
arr_x = np.array(donnees_x, dtype=float)
u_x = np.array(incertitudes_x, dtype=float)
arr_y = np.array(donnees_y, dtype=float)
u_y = np.array(incertitudes_y, dtype=float)

# Premier ajustement linéaire classique (sans poids) pour estimer la pente brute
pente_initiale, _ = np.polyfit(arr_x, arr_y, 1)

# Calcul du facteur de pondération statistique combinant les incertitudes en X et en Y
poids_stat = 1.0 / np.sqrt(u_y**2 + (pente_initiale * u_x) ** 2)

# Ajustement linéaire définitif par la méthode des moindres carrés pondérés (Chi-deux)
pente_khi2, ordonnee_khi2 = np.polyfit(arr_x, arr_y, 1, w=poids_stat)

# Initialisation du moteur de génération de nombres aléatoires de NumPy
generateur_rnd = np.random.default_rng()

# Création des matrices de tirages aléatoires (Lois normales centrées sur chaque point)
simuls_y = generateur_rnd.normal(arr_y, u_y, (TOTAL_TIRAGES, len(arr_y)))
simuls_x = generateur_rnd.normal(arr_x, u_x, (TOTAL_TIRAGES, len(arr_x)))

# Initialisation de listes vides pour stocker le résultat de chaque simulation
pentes_calculees = []
ordonnees_calculees = []

# --- BOUCLE DE SIMULATION MONTE-CARLO ---
# On reproduit fidèlement la structure itérative ligne par ligne
for i in range(TOTAL_TIRAGES):
    # Évaluation de la pondération spécifique à chaque set de données simulé
    poids_local = 1.0 / np.sqrt(u_y**2 + (pente_khi2 * u_x) ** 2)

    # Ajustement de la droite sur les points perturbés du tirage en cours
    coefficients_i = np.polyfit(simuls_x[i], simuls_y[i], 1, w=poids_local)

    # Sauvegarde de la pente (index 0) et de l'ordonnée à l'origine (index 1) obtenues
    pentes_calculees.append(coefficients_i[0])
    ordonnees_calculees.append(coefficients_i[1])

# Conversion des listes de résultats en tableaux NumPy pour l'analyse statistique
pentes_calculees = np.array(pentes_calculees)
ordonnees_calculees = np.array(ordonnees_calculees)

# Extraction des écarts-types des distributions (Incertitudes-types de Monte-Carlo)
ecart_type_a = np.std(pentes_calculees)
ecart_type_b = np.std(ordonnees_calculees)

# Calcul des incertitudes élargies avec la statistique de Student
U_pente = FACTEUR_STUDENT * ecart_type_a
U_ordonnee = FACTEUR_STUDENT * ecart_type_b

# Évaluation de la pertinence du modèle (Calcul du coefficient R²)
y_estimations = pente_khi2 * arr_x + ordonnee_khi2
somme_carres_residus = np.sum((arr_y - y_estimations) ** 2)
somme_carres_totaux = np.sum((arr_y - np.mean(arr_y)) ** 2)
coefficient_R2 = 1 - (somme_carres_residus / somme_carres_totaux)

# =============================================================================
# RESTITUTION DES RÉSULTATS DANS LA CONSOLE
# =============================================================================

print("~" * 60)
print(f" COMPTE-RENDU STATISTIQUE : {TITRE_MANIP} ")
print("~" * 60)
print(f" Formule de la droite : Y = a * X + b")
print(f" Pente calculée (a)   : {pente_khi2:.4E}")
print(f" Incertitude U(a)     : {U_pente:.2E} (Confiance 95 %)")
print(f" Ordonnée à l'orig.(b): {ordonnee_khi2:.4E}")
print(f" Incertitude U(b)     : {U_ordonnee:.2E} (Confiance 95 %)")
print(f" Coefficient R²       : {coefficient_R2:.5f}")
print(f" Itérations MC        : {TOTAL_TIRAGES}")
print("~" * 60)

# =============================================================================
# GÉNÉRATION DES PANELS GRAPHIQUES
# =============================================================================

# Création d'un espace graphique composé de 3 sous-figures alignées horizontalement
fig, zones = plt.subplots(1, 3, figsize=(15, 4.8))

# --- Graphe 1 : Droite d'ajustement et données réelles avec barres d'erreur ---
zones[0].errorbar(
    arr_x,
    arr_y,
    xerr=u_x,
    yerr=u_y,
    fmt="o",
    color="#0066cc",
    ecolor="#666666",
    capsize=3,
    label="Mesures avec incertitudes",
)
x_espace = np.linspace(arr_x.min() * 0.9, arr_x.max() * 1.1, 150)
zones[0].plot(
    x_espace,
    pente_khi2 * x_espace + ordonnee_khi2,
    color="#cc0000",
    linewidth=1.8,
    label=f"Modèle régressé",
)
zones[0].set_xlabel(NOM_AXE_X)
zones[0].set_ylabel(NOM_AXE_Y)
zones[0].set_title("Régression Linéaire Pondérée")
zones[0].legend(fontsize=9)
zones[0].grid(True, linestyle=":", alpha=0.5)

# --- Graphe 2 : Distribution de la pente (a) par Monte-Carlo ---
zones[1].hist(pentes_calculees, bins=70, color="#0066cc", alpha=0.7, density=True)
zones[1].axvline(pente_khi2, color="#cc0000", linestyle="-", label=f"a = {pente_khi2:.3E}")
zones[1].axvline(
    pente_khi2 - U_pente,
    color="#ff9900",
    linestyle="--",
    label=f"± U(a) = {U_pente:.2E}",
)
zones[1].axvline(pente_khi2 + U_pente, color="#ff9900", linestyle="--")
zones[1].set_xlabel("Valeurs de la pente a")
zones[1].set_title("Distribution Monte-Carlo de a")
zones[1].legend(fontsize=9)
zones[1].grid(True, linestyle=":", alpha=0.5)

# --- Graphe 3 : Distribution de l'ordonnée à l'origine (b) par Monte-Carlo ---
zones[2].hist(
    ordonnees_calculees, bins=70, color="#2ca02c", alpha=0.7, density=True
)
zones[2].axvline(ordonnee_khi2, color="#cc0000", linestyle="-", label=f"b = {ordonnee_khi2:.3E}")
zones[2].axvline(
    ordonnee_khi2 - U_ordonnee,
    color="#ff9900",
    linestyle="--",
    label=f"± U(b) = {U_ordonnee:.2E}",
)
zones[2].axvline(ordonnee_khi2 + U_ordonnee, color="#ff9900", linestyle="--")
zones[2].set_xlabel("Valeurs de l'ordonnée b")
zones[2].set_title("Distribution Monte-Carlo de b")
zones[2].legend(fontsize=9)
zones[2].grid(True, linestyle=":", alpha=0.5)

# Ajustement global de la géométrie des graphiques pour éviter les chevauchements
plt.tight_layout()
plt.show()
 