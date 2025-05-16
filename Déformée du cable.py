# ===============================
#    MODELE DE DEFORMEE DE CABLE
# ===============================

# --- 1. IMPORTS ---
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

# --- 2. CONSTANTES PHYSIQUES ET PARAMÈTRES ---
#region ---Grandeurs physiques---
g = 9.81
rho_eau = 1025  # kg/m³
mu = 1.1*10**(-3)
#endregion

#region --- CÂBLE ---
L_tot = 50.0       # Longueur totale du câble (m)
n = 150            # Nombre de segments
delta_s = L_tot / n# Longueur d'un segement

D = 0.0063              # Diamètre du câble (m)
S = (D**2) * np.pi / 4  # section du cable
mass_lin = 0.035        # Masse linéique du câble (kg/m)

P_c = - (mass_lin - rho_eau * (np.pi * D ** 2 / 4)) * g * delta_s  # Poids apparent
#endregion

#region --- bateau ---
v_kt = 3 # Vitesse en noeud
v_bateau = v_kt * 0.5144  # Vitesse en m/s
#endregion

#region --- SONAR ---
m_sondeur_imm = 6.7             # Masse apparente (plongée)
D_sondeur = 0.06
S_sondeur = (D_sondeur**2) * np.pi / 4  # Surface frontale
L_sondeur = 0.85
L_ail = 0.105
h_ail = 0.028

S_mouillee = np.pi * D_sondeur * (L_sondeur + D_sondeur) + 6 * L_ail * h_ail

#endregion

#region --- COEFFICIENTS HYDRODYNAMIQUES ---



Cd_sondeur = 1.0 # Trainée tangentielle sondeur
Re = rho_eau * v_bateau * L_sondeur/mu
Cf_sondeur = 0.455/((np.log(Re))**2.58) # trainée dynamique de frottement visqueux

Cd = 1.01     # Traînée normale section cable

#endregion

#region --- 4. FORCES ----

#force sur l_element
def calcul_forces(theta, L_tot = L_tot, v_kt = v_kt, n = n):

    v_bateau = v_kt * 0.5144  # Vitesse en m/s
    delta_s = L_tot / n
    v_n = v_bateau * np.cos(theta)
    v_t = v_bateau * np.sin(theta)

    Re = rho_eau * abs(v_t) * delta_s / mu
    Cf = 0.455 / ((np.log(Re)) ** 2.58)

    F_n = - 0.5 * rho_eau * v_n**2 * Cd * D * delta_s
    F_t = - 0.5 * rho_eau * v_t**2 * Cf * D * delta_s * np.pi

    return F_n, F_t

print(f"Pour un angle faible (10^-4) l'effort sur le cable est \nFn = {calcul_forces(0.0001)[0]} N et Ft =  {calcul_forces(0.0001)[1]} N."
      f"\nPour l'angle maximum immaginable de l'effort du sondeur sur le premier élement de cable "
      f"(la trainée s'opposerai à l'inclinaison d'où angle max) "
      f"on aurait \nFn = {calcul_forces(10)[0]} N et Ft = {calcul_forces(10)[1]} N"
      f"\nOn peud encadrer au moins par ces angles l'inclinaison du premier element de câble")

def effort_sondeur(v = v_kt):
    v_bateau = v * 0.5144
    Re = rho_eau * v_bateau * L_sondeur / mu
    Cf_sondeur = 0.455 / ((np.log(Re)) ** 2.58)  # trainée dynamique de frottement visqueux

    Fn = m_sondeur_imm * g
    Ft = 0.5 * rho_eau * (Cd_sondeur * S_sondeur + Cf_sondeur * S_mouillee) * v_bateau**2
    return (Ft, -Fn)


#endregion

#region --- 5. CALCUL DE LA DÉFORMÉE ---

def resoudre_deformee(L_tot = L_tot, v_kt = v_kt, n = n) :

    delta_s = L_tot / n

    alpha = np.zeros(n+1)
    #Beta = alpha[i+1]
    x = np.zeros(n + 1)
    z = np.zeros(n + 1)
    F = []

    eff_sond = effort_sondeur(v_kt)
    F.append(np.sqrt(eff_sond[0] ** 2 + eff_sond[1] ** 2))
    x[0] = 0
    z[0] = 0


    ## on initialise l'angle de la première section
    alpha[0] = np.arctan2(eff_sond[0],eff_sond[1])

    F_nc, F_tc = calcul_forces(alpha[0], L_tot = L_tot, v_kt = v_kt, n = n)

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

        F_nc, F_tc = calcul_forces(alpha[i], L_tot, v_kt, n)

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
#endregion

#region --- 6. AFFICHAGE DEFORMEE ---

def tracer_deformee(x, z):
    plt.figure(figsize=(10, 6))
    plt.plot(x, z)  # z orienté vers le bas
    plt.xlabel("Distance horizontale (m)")
    plt.ylabel("Profondeur (m)")
    plt.title("Déformée du câble")
    plt.grid(True)
    plt.axis("equal")
    plt.show()
#endregion

# --- 7. MAIN ---
if __name__ == "__main__":
    x, z = resoudre_deformee()
    tracer_deformee(x, z)

    layback = abs(x[0])   # distance horizontale depuis le bateau (en 0)
    profondeur = abs(z[0])  # profondeur positive

    print(f"\n--- RÉSULTATS ---")
    print(f"Layback (distance horizontale sonar-bateau) : {layback:.2f} m")
    print(f"Profondeur du sonar : {profondeur:.2f} m")


    vitesses = [0.5, 1, 1.5, 2, 3, 5]
    longueurs = [10, 25, 50, 100, 150, 200]

    laybacks = np.zeros((len(vitesses), len(longueurs)))
    profondeurs = np.zeros_like(laybacks)

    profondeur_constructeur = np.array([
        [10, 9, 9, 9, 8, 5],
        [24, 24, 22, 20, 14, 8],
        [49, 46, 39, 31, 20, 10],
        [98, 81, 58, 43, 25, 12],
        [143, 105, 71, 50, 28, 13],
        [186, 124, 80, 55, 30, 14]
    ])

    layback_constructeur = np.array([
        [0, 0, 1, 2, 4, 7],
        [1, 4, 8, 12, 18, 22],
        [4, 15, 27, 35, 42, 47],
        [17, 51, 73, 83, 92, 97],
        [37, 95, 121, 133, 142, 147],
        [63, 141, 171, 183, 192, 197]
    ])

    for i, v in enumerate(vitesses):
        for j, L in enumerate(longueurs):
            x, z = resoudre_deformee(L_tot=L, v_kt=v)
            laybacks[i, j] = abs(x[0])
            profondeurs[i, j] = abs(z[0])

    ecarts_layback = (laybacks - layback_constructeur)
    ecarts_profondeur = (profondeurs - profondeur_constructeur)

    df_prof = pd.DataFrame(profondeur_constructeur, index=longueurs, columns=vitesses)
    df_layb = pd.DataFrame(layback_constructeur, index=longueurs, columns=vitesses)

    plt.figure(figsize=(10, 6))
    for vitesse in df_prof.columns:
        plt.plot(df_prof.index, df_prof[vitesse], marker='o', label=f"{vitesse} kt")

    plt.title("Profondeur en fonction de la longueur de câble")
    plt.xlabel("Longueur de câble (m)")
    plt.ylabel("Profondeur (m)")
    plt.grid(True)
    plt.legend(title="Vitesse")
    plt.tight_layout()
    plt.show()

    plt.figure(figsize=(10, 6))
    for vitesse in df_layb.columns:
        plt.plot(df_layb.index, df_layb[vitesse], marker='s', label=f"{vitesse} kt")

    plt.title("Layback en fonction de la longueur de câble")
    plt.xlabel("Longueur de câble (m)")
    plt.ylabel("Layback (m)")
    plt.grid(True)
    plt.legend(title="Vitesse")
    plt.tight_layout()
    plt.show()


    # --- ÉCARTS RELATIFS ---
    laybacks_rel = laybacks / np.array(longueurs)[:, None]
    profondeurs_rel = profondeurs / np.array(longueurs)[:, None]
    laybacks_ref_rel = layback_constructeur / np.array(longueurs)[:, None]
    profondeurs_ref_rel = profondeur_constructeur / np.array(longueurs)[:, None]

    ecarts_rel_layback = laybacks_rel - laybacks_ref_rel
    ecarts_rel_profondeur = profondeurs_rel - profondeurs_ref_rel

    # --- AFFICHAGE TABLEAU ---

    print("\nLAYBACK RELATIFS DE REF (Layback / Longueur câble) :")
    print(laybacks_ref_rel.round(3))

    print("\nPROFONDEURS RELATIVES DE REF (Profondeur / Longueur câble) :")
    print(profondeurs_ref_rel.round(3))

    print("\nLAYBACK RELATIFS (Layback / Longueur câble) :")
    print(laybacks_rel.round(3))

    print("\nPROFONDEURS RELATIVES (Profondeur / Longueur câble) :")
    print(profondeurs_rel.round(3))

    df_rel_layback = pd.DataFrame(ecarts_rel_layback, index=longueurs, columns=vitesses)
    df_rel_prof = pd.DataFrame(ecarts_rel_profondeur, index=longueurs, columns=vitesses)

    print("\nÉCARTS RELATIFS (Layback / Longueur câble) :")
    print(df_rel_layback.round(3))

    print("\nÉCARTS RELATIFS (Profondeur / Longueur câble) :")
    print(df_rel_prof.round(3))

    # --- TRACÉ ---
    plt.figure(figsize=(10, 6))
    for v in vitesses:
        plt.plot(longueurs, df_rel_layback[v], marker='o', label=f"{v} kt")
    plt.title("Écart relatif de layback (modèle - constructeur)")
    plt.xlabel("Longueur de câble (m)")
    plt.ylabel("Écart relatif (layback / L)")
    plt.grid(True)
    plt.legend(title="Vitesse")
    plt.tight_layout()
    plt.show()

    plt.figure(figsize=(10, 6))
    for v in vitesses:
        plt.plot(longueurs, df_rel_prof[v], marker='s', label=f"{v} kt")
    plt.title("Écart relatif de profondeur (modèle - constructeur)")
    plt.xlabel("Longueur de câble (m)")
    plt.ylabel("Écart relatif (profondeur / L)")
    plt.grid(True)
    plt.legend(title="Vitesse")
    plt.tight_layout()
    plt.show()
