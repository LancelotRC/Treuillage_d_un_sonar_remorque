# ===============================
#    MODELE DE DEFORMEE DE CABLE
# ===============================
from urllib.request import FTPHandler

# --- 1. IMPORTS ---
import numpy as np
import matplotlib.pyplot as plt

# --- 2. CONSTANTES PHYSIQUES ET PARAMÈTRES ---
g = 9.81
rho_eau = 1025  # kg/m³

# --- CÂBLE ---
L_tot = 50.0       # Longueur totale du câble (m)
n = 150            # Nombre de segments
delta_s = L_tot / n# Longueur d'un segement

D = 0.0063              # Diamètre du câble (m)
S = (D**2) * np.pi / 4  # section du cable
mass_lin = 0.          # Masse linéique du câble (kg/m)

P_c = - (mass_lin - rho_eau * (np.pi * D ** 2 / 4)) * g * delta_s  # Poids apparent

# --- bateau ---
v_kt = 3                  # Vitesse en noeud
v_bateau = v_kt * 0.5144  # Vitesse en m/s



# --- SONAR ---
m_sondeur_imm = 6.7             # Masse apparente (plongée)
D_sondeur = 0.06
S_sondeur = (D_sondeur**2) * np.pi / 4  # Surface frontale


# --- COEFFICIENTS HYDRODYNAMIQUES ---
Cd_sondeur = 0.7 # Trainée tangentielle sondeur

Cd = 1.01     # Traînée normale section cable
Cf = 0.01    # Traînée tangentielle section cable

# --- 4. FORCES ----

#force sur l_element
def calcul_forces(theta):
    delta_s = L_tot / n
    v_n = v_bateau * np.cos(theta)
    v_t = v_bateau * np.sin(theta)

    F_n = - 0.5 * rho_eau * v_n**2 * Cd * D * delta_s
    F_t = - 0.5 * rho_eau * v_t**2 * Cf * D * delta_s * np.pi

    return F_n, F_t

def effort_sondeur():
    Fn = m_sondeur_imm * g
    Ft = 0.5 * rho_eau * Cd_sondeur * S_sondeur * v_bateau**2
    return (-Ft, -Fn)


# --- 5. CALCUL DE LA DÉFORMÉE ---

def resoudre_deformee():
    alpha = np.zeros(n+1)
    #Beta = alpha[i+1]
    x = np.zeros(n + 1)
    z = np.zeros(n + 1)
    F = []

    eff_sond = effort_sondeur()
    F.append(np.sqrt(eff_sond[0] ** 2 + eff_sond[1] ** 2))
    x[0] = 0
    z[0] = 0


    ## on initialise l'angle de la première section
    alpha[0] = np.arctan2(eff_sond[0],eff_sond[1])

    F_nc, F_tc = calcul_forces(alpha[0])

    #region == calcul ==
        #projection sur la normal du câble :
    # -F_nc + np.sin(alpha[1]) * F[1] - np.sin(alpha[0]) * (F[0] - P_c) = 0

        #projection sur la tengente du câble :
    # -F_tc + np.cos(alpha[1]) * F[1] - np.cos(alpha[0]) * (F[0] + P_c) = 0

    # (np.sin(alpha[1]) * F[1])²  = (F_nc + np.sin(alpha[0]) * (F[0] - P_c))²
    # (np.cos(alpha[1]) * F[1])²  = (F_nc + np.cos(alpha[0]) * (F[0] + P_c))²

    # ( (np.sin(alpha[1])² + (np.cos(alpha[1])² ) x F[1]²
    #        = (F_nc + np.sin(alpha[0]) * (F[0] - P_c))² + (F_nc + np.cos(alpha[0]) * (F[0] + P_c))²

    # F[1] = rac ( [ (F_nc + np.sin(alpha[0]) * (F[0] - P_c))² + (F_nc + np.cos(alpha[0]) * (F[0] + P_c))² ] )
    # alpha[1] = arcsin(-(1 / F[1]) * (F_nc + sin(-alpha[0]) * (F[0] + P_c)))
    #endregion

    F_i = np.sqrt(
        (F_nc - np.sin(alpha[0]) * (F[0] - P_c)) ** 2 +
        (F_tc - np.cos(alpha[0]) * (F[0] - P_c)) ** 2
    )

    alpha[1] = np.arctan2(
        -(F_nc - np.sin(alpha[0]) * (F[0] - P_c)),
        -(F_tc - np.cos(alpha[0]) * (F[0] - P_c))
    )
    F.append(F_i)


    for i in range(n):
        x[i+1] = x[i] + delta_s * np.sin(alpha[i])
        z[i+1] = z[i] - delta_s * np.cos(alpha[i])

        F_nc, F_tc = calcul_forces(alpha[i])

        F_imoins1 = F[i]

        F_i = np.sqrt(
            (F_nc - np.sin(alpha[i]) * (F_imoins1 - P_c)) ** 2 +
            (F_tc - np.cos(alpha[i]) * (F_imoins1 - P_c)) ** 2
        )

        # on actualise l'angle et on ajoute la positionn
        alpha[i + 1] = np.arctan2(
            -(F_nc - np.sin(alpha[i]) * (F_imoins1 - P_c)),
            -(F_tc - np.cos(alpha[i]) * (F_imoins1 - P_c))
        )
        F.append(F_i)
    return x - x[-1], z - z[-1]


# --- 6. AFFICHAGE ---

def tracer_deformee(x, z):
    plt.figure(figsize=(10, 6))
    plt.plot(x, z)  # z orienté vers le bas
    plt.xlabel("Distance horizontale (m)")
    plt.ylabel("Profondeur (m)")
    plt.title("Déformée du câble")
    plt.grid(True)
    plt.axis("equal")
    plt.show()


# --- 7. MAIN ---
if __name__ == "__main__":
    x, z = resoudre_deformee()
    tracer_deformee(x, z)

    layback = abs(x[0])   # distance horizontale depuis le bateau (en 0)
    profondeur = abs(z[0])  # profondeur positive

    print(f"\n--- RÉSULTATS ---")
    print(f"Layback (distance horizontale sonar-bateau) : {layback:.2f} m")
    print(f"Profondeur du sonar : {profondeur:.2f} m")
