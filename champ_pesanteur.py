import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider

# --- 1. Données Physiques Réelles de la Terre (SI) ---
G = 6.6743e-11
M_terre = 5.972e24
R_terre = 6.371e6      # Rayon moyen (m)
omega_terre = 7.292e-5 # Vitesse de rotation réelle (rad/s)

# --- 2. Calculs Physiques Exacts ---
g_grav_norme = (G * M_terre) / (R_terre**2)
epsilon_reel = (omega_terre**2 * R_terre) / g_grav_norme

def calcul_champs(latitude_deg):
    """
    Calcule les vecteurs g_grav, a_centrifuge et g_pesanteur réels
    et applique un facteur d'exagération VISUELLE pour la géométrie de la loupe.
    """
    lat_rad = np.radians(latitude_deg)
    
    # 1. Vecteur gravitationnel pur (vers le centre de la Terre -> -X)
    g_grav = np.array([-g_grav_norme, 0.0])
    
    # 2. Accélération d'entraînement centrifuge réelle
    H = R_terre * np.cos(lat_rad)
    a_e_norme_reelle = (omega_terre**2) * H
    a_e_reelle = np.array([a_e_norme_reelle * np.cos(lat_rad), -a_e_norme_reelle * np.sin(lat_rad)])
    
    # 3. Champ de pesanteur réel résultant
    g_pes_reelle = g_grav + a_e_reelle
    g_pes_norme = np.linalg.norm(g_pes_reelle)
    
    # 4. Angle de déviation réel de la verticale (en degrés)
    angle_rad = np.abs(np.arctan2(g_pes_reelle[1], g_pes_reelle[0])) - np.pi
    angle_deg_reel = np.degrees(np.abs(angle_rad + np.pi))
    if latitude_deg == 0 or latitude_deg == 90:
        angle_deg_reel = 0.0
        
    # --- EXAGÉRATION VISUELLE (Facteur 30) POUR LA LOUPE ---
    facteur_zoom = 30.0
    a_e_visuelle = a_e_reelle * facteur_zoom
    g_pes_visuelle = g_grav + a_e_visuelle
    
    return g_grav, a_e_visuelle, g_pes_visuelle, g_pes_norme, angle_deg_reel

# --- 3. Pré-calcul de la courbe du haut (Fonction de la latitude) ---
latitudes_axes = np.linspace(0, 90, 200)
g_pes_liste = []

for lat in latitudes_axes:
    _, _, _, g_p, _ = calcul_champs(lat)
    g_pes_liste.append(g_p)

# --- 4. Initialisation de la Figure (Structure 2x2 conservée) ---
fig = plt.figure(figsize=(14, 8))
plt.subplots_adjust(bottom=0.22, wspace=0.3, hspace=0.35)

ax_globe = plt.subplot(2, 2, 1)
ax_zoom = plt.subplot(2, 2, 3)
ax_globe.set_aspect('equal')
ax_zoom.set_aspect('equal')

ax_g_curve = plt.subplot(2, 2, 2)
ax_vide = plt.subplot(2, 2, 4) 
ax_vide.axis('off')           

# Tracé de la courbe en haut à droite
line_g_curve, = ax_g_curve.plot(latitudes_axes, g_pes_liste, color='crimson', lw=2)
marker_g_curve, = ax_g_curve.plot([], [], 'ro', markersize=7)

# --- 5. Dessin du Globe Statique ---
theta = np.linspace(0, 2*np.pi, 100)
ax_globe.plot(R_terre * np.cos(theta), R_terre * np.sin(theta), color='lightblue', lw=1.5)
ax_globe.fill(R_terre * np.cos(theta), R_terre * np.sin(theta), color='lightblue', alpha=0.3)
ax_globe.axvline(0, color='blue', linestyle='--', alpha=0.5)
ax_globe.axhline(0, color='grey', linestyle=':', alpha=0.5)

point_globe, = ax_globe.plot([], [], 'go', markersize=6)
loupe_cercle = plt.Circle((0, 0), R_terre*0.2, color='black', fill=False, linestyle='--', alpha=0.5)
ax_globe.add_patch(loupe_cercle)

# Éléments de la Loupe
point_zoom, = ax_zoom.plot(0, 0, 'go', markersize=8, zorder=5)

quiver_grav = ax_zoom.quiver(0, 0, 0, 0, angles='xy', scale_units='xy', scale=1, color='teal', label=r'$\vec{g}_{grav}$')
quiver_cent = ax_zoom.quiver(0, 0, 0, 0, angles='xy', scale_units='xy', scale=1, color='gold', label=r'$30 \times \vec{a}_e$ (exagéré)')
quiver_pes = ax_zoom.quiver(0, 0, 0, 0, angles='xy', scale_units='xy', scale=1, color='crimson', label=r'$\vec{g}_{visuel}$')

# Texte brut écrit en VERT SAPIN (forestgreen) et en gras, placé en haut du cadre
text_epsilon = ax_zoom.text(-11.5, 1.2, "", fontsize=10, color='forestgreen', fontweight='bold')

# --- 6. Fonction de mise à jour ---
def update(val):
    lat_deg = slider_lat.val
    lat_rad = np.radians(lat_deg)
    
    g_grav, a_e_vis, g_pes_vis, g_p_n, ang_d = calcul_champs(lat_deg)
    
    # 1. Point sur le globe complet
    x_g = R_terre * np.cos(lat_rad)
    y_g = R_terre * np.sin(lat_rad)
    point_globe.set_data([x_g], [y_g])
    loupe_cercle.set_center((x_g, y_g))
    
    # 2. Flèches dans la loupe
    global quiver_grav, quiver_cent, quiver_pes
    quiver_grav.set_UVC(g_grav[0], g_grav[1])
    quiver_cent.set_UVC(a_e_vis[0], a_e_vis[1])
    quiver_pes.set_UVC(g_pes_vis[0], g_pes_vis[1])
    
    # Texte mis à jour en vert
    text_epsilon.set_text(f"Ordre de grandeur :  \u03b5 = {epsilon_reel:.5f}\n"
                          f"Déviation réelle :  \u03b1 = {ang_d:.3f}°\n"
                          f"Intensité locale :  g = {g_p_n:.3f} m/s\u00b2")
    
    # 3. Marqueur sur la courbe
    marker_g_curve.set_data([lat_deg], [g_p_n])
    
    fig.canvas.draw_idle()

# --- 7. Création du Curseur ---
ax_slider = plt.axes([0.15, 0.06, 0.7, 0.03])
slider_lat = Slider(ax=ax_slider, label='Latitude (\u03bb) ', valmin=0.0, valmax=90.0, valinit=45.0, valfmt='%1.1f°', color='seagreen')
slider_lat.on_changed(update)

# --- 8. Habillage ---
ax_globe.set_xlim(-R_terre*1.2, R_terre*1.2)
ax_globe.set_ylim(-R_terre*1.2, R_terre*1.2)
ax_globe.set_title("Position sur la Terre", fontsize=11, fontweight='bold')
ax_globe.axis('off')

# Échelle fixe de la loupe
ax_zoom.set_xlim(-12, 2)
ax_zoom.set_ylim(-3, 3)
ax_zoom.set_title("Loupe (Forces en $m\\cdot s^{-2}$ - Effet centrifuge $\\times 30$)", fontsize=11, fontweight='bold')
ax_zoom.grid(True, linestyle=':', alpha=0.6)
ax_zoom.legend(loc='lower left', fontsize=9)
ax_zoom.set_xlabel("Axe Radial (vers l'extérieur)")
ax_zoom.set_ylabel("Axe Tangentiel (vers le Nord)")

ax_g_curve.set_ylabel(r"Pesanteur $g$ ($m\cdot s^{-2}$)")
ax_g_curve.set_title("Intensité du champ de pesanteur $g(\\lambda)$", fontsize=11, fontweight='bold')
ax_g_curve.grid(True, linestyle='--', alpha=0.5)
ax_g_curve.set_xlabel("Latitude \u03bb (degrés)")

update(None)
plt.show()