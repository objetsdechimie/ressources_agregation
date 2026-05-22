# -*- coding: utf-8 -*-
"""
Densités de probabilité radiale des orbitales atomiques 1s, 2s, 3s
r est adimensionné par le rayon de Bohr a0 = 53 pm
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import quad

# ─────────────────────────────────────────────
# PARTIES RADIALES R(r) ET DENSITÉS D(r) = r²·R²
# ─────────────────────────────────────────────

ORBITALES = {
    (1, 0): {"nom": "1s", "couleur": "#2ecc71", "style": "-."},
    (2, 0): {"nom": "2s", "couleur": "#3498db", "style": "--"},
    (3, 0): {"nom": "3s", "couleur": "#e74c3c", "style": "-"},
}

def partie_radiale(n: int, l: int, r: np.ndarray) -> np.ndarray:
    """Retourne la partie radiale R_{n,l}(r) normalisée."""
    if n == 1 and l == 0:   # 1s
        return 2 * np.exp(-r)
    if n == 2 and l == 0:   # 2s  — CORRECTION : 1/(2√2) · (2 - r) · e^(-r/2)
        return (1 / (2 * np.sqrt(2))) * (2 - r) * np.exp(-r / 2)
    if n == 3 and l == 0:   # 3s  — CORRECTION : polynôme exact
        return (2 / (81 * np.sqrt(3))) * (27 - 18*r + 2*r**2) * np.exp(-r / 3)
    raise ValueError(f"Orbitale (n={n}, l={l}) non définie.")

def densite_radiale(n: int, l: int, r: np.ndarray) -> np.ndarray:
    """Retourne la densité de probabilité radiale D(r) = r² · R²(r)."""
    R = partie_radiale(n, l, r)
    return r**2 * R**2

# ─────────────────────────────────────────────
# RAYON LE PLUS PROBABLE  (maximum de D(r))
# ─────────────────────────────────────────────

def rayon_max(n: int, l: int, r: np.ndarray) -> float:
    """Retourne le rayon correspondant au maximum de la densité radiale."""
    D = densite_radiale(n, l, r)
    return r[np.argmax(D)]   # np.argmax : plus robuste que la boucle while

# ─────────────────────────────────────────────
# RAYON MOYEN  <r> = ∫ r · D(r) dr
# ─────────────────────────────────────────────

def rayon_moyen(n: int, l: int) -> float:
    """Calcule <r> = ∫₀^∞ r · D(r) dr par intégration numérique."""
    integrand = lambda r: r * densite_radiale(n, l, r)
    valeur, _ = quad(integrand, 0, np.inf)
    return valeur

# ─────────────────────────────────────────────
# VÉRIFICATION DE LA NORMALISATION ∫ D(r) dr = 1
# ─────────────────────────────────────────────

def verifier_normalisation(n: int, l: int) -> float:
    """Retourne ∫₀^∞ D(r) dr (doit être ≈ 1.0)."""
    integrand = lambda r: densite_radiale(n, l, r)
    valeur, _ = quad(integrand, 0, np.inf)
    return valeur

# ─────────────────────────────────────────────
# TRACÉ
# ─────────────────────────────────────────────

def tracer_densites(orbitales: dict, r: np.ndarray) -> None:
    fig, ax = plt.subplots(figsize=(9, 5))

    for (n, l), props in orbitales.items():
        D = densite_radiale(n, l, r)
        ax.plot(r, D,
                color=props["couleur"],
                linestyle=props["style"],
                linewidth=2,
                label=f'OA {props["nom"]}')

        # Marqueur du maximum
        r_pic = rayon_max(n, l, r)
        D_pic = densite_radiale(n, l, r_pic)
        ax.axvline(r_pic, color=props["couleur"], linewidth=0.8, alpha=0.4)
        ax.scatter([r_pic], [D_pic], color=props["couleur"], zorder=5, s=50)

    ax.set_xlabel(r"$r \ / \ a_0$", fontsize=13)
    ax.set_ylabel(r"$D(r) = r^2 \, R_{nl}^2(r)$", fontsize=13)
    ax.set_title("Densité de probabilité radiale — orbitales ns (H)", fontsize=14)
    ax.legend(fontsize=11)
    ax.set_xlim(0, r[-1])
    ax.set_ylim(bottom=0)
    ax.grid(True, linestyle=":", alpha=0.5)
    fig.tight_layout()
    plt.savefig("/mnt/user-data/outputs/densites_radiales.png", dpi=150)
    plt.show()

# ─────────────────────────────────────────────
# RÉSULTATS NUMÉRIQUES
# ─────────────────────────────────────────────

def afficher_resultats(orbitales: dict, r: np.ndarray) -> None:
    print(f"\n{'OA':<6} {'r_max (a0)':>12} {'<r> (a0)':>12} {'Norm.':>10}")
    print("─" * 44)
    for (n, l), props in orbitales.items():
        rmax  = rayon_max(n, l, r)
        rmoy  = rayon_moyen(n, l)
        norme = verifier_normalisation(n, l)
        print(f"{props['nom']:<6} {rmax:>12.2f} {rmoy:>12.2f} {norme:>10.4f}")
    print()

# ─────────────────────────────────────────────
# MAIN
# ─────────────────────────────────────────────

if __name__ == "__main__":
    r = np.linspace(0, 35, 2000)   # étendu à 35 a0 pour bien voir la 3s
    tracer_densites(ORBITALES, r)
    afficher_resultats(ORBITALES, r)

