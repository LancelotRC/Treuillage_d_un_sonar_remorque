# =======================================
# ===   MODELE DE DEFORMEE DE CABLE   ===
# =======================================

#region --- 1. IMPORTS ---
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
#endregion

#region --- 2. CONSTANTES PHYSIQUES ET PARAMÈTRES ---
#region ---Grandeurs physiques---
g = 9.81
rho_eau = 1025  # kg/m³
mu = 1.1*10**(-3)
#endregion

#region --- CÂBLE ---
L_tot = 50.0            # Longueur totale du câble (m)
n = 5000                 # Nombre de segments
D = 0.0063              # Diamètre du câble (m)
mass_lin = 0.035        # Masse linéique du câble (kg/m)

delta_s = L_tot / n     # Longueur d'un segement
S = (D**2) * np.pi / 4  # section du cable
P_c = - (mass_lin - rho_eau * (np.pi * D ** 2 / 4)) * g * delta_s
                        # Poids apparent
#endregion

#region --- BATEAU ---
v_kt = 3 # Vitesse en noeud
v_bateau = v_kt * 0.5144  # Vitesse en m/s
#endregion

#region --- SONAR ---
m_sondeur_imm = 6.7     # Masse immergée
D_sondeur = 0.06        # Diametre
L_sondeur = 0.85        # Longueur
L_ail = 0.105           # longueur des ailettes
h_ail = 0.028           # hauteur radiale des ailettes

S_sondeur = (D_sondeur**2) * np.pi / 4
                        # Surface frontale
S_mouillee = np.pi * D_sondeur * (L_sondeur + D_sondeur) + 6 * L_ail * h_ail
                        # Surface au contact de l'eau

#endregion

#region --- COEFFICIENTS HYDRODYNAMIQUES ---

Cd_sondeur = 1.0        # Trainée tangentielle sondeur
Cd = 1.01               # Traînée normale section cable

Re = rho_eau * v_bateau * L_sondeur/mu
                        # Reynolds pour le sondeur
Cf_sondeur = 0.455/((np.log(Re))**2.58)
                        # trainée dynamique de frottement visqueux
#endregion
#endregion

#region --- 4. FORCES ----

def calcul_forces(theta, L_tot = L_tot, v_kt = v_kt, n = n):
    # calcul des efforts sur une section du cable
    v_bateau = v_kt * 0.5144  # Vitesse en m/s
    dl = L_tot / n
    v_n = v_bateau * np.cos(theta)
    v_t = v_bateau * np.sin(theta)

    Re = rho_eau * abs(v_t) * dl / mu
    Cf = 0.455 / ((np.log(Re)) ** 2.58)

    F_n = - 0.5 * rho_eau * v_n**2 * Cd * D * dl
    F_t = - 0.5 * rho_eau * v_t**2 * Cf * D * dl * np.pi

    return F_n, F_t

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

    dl = L_tot / n

    #Beta = alpha[i+1]
    alpha = np.zeros(n+1)
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
        x[i+1] = x[i] + dl * np.sin(alpha[i])
        z[i+1] = z[i] - dl * np.cos(alpha[i])

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

    # --- PARAMÈTRES DE CONTRÔLE ---
    afficher_encadrement_angles_et_efforts = False  # affiche le print ci-dessous qui justifie le choix de l'angle initial
    afficher_deformee = True                        # affiche la déformée calculée pour une valeur (réglage dans la Partie 1 du main)
    afficher_graphiques_comparaison = False         # Permet de montrer sur un graph les résultats (pas très parlant)
    afficher_tableaux_ecarts = True                 # Renvoie les tableaux numériques des écarts absolus et relatifs (Partie 4)
    afficher_graph = True                           # Trace le graph des écarts relatifs au constructeur pour profondeur et Layback

    if afficher_encadrement_angles_et_efforts:
        print(
            f"Pour un angle faible (10^-4) l'effort sur le cable est \nFn = {calcul_forces(0.0001)[0]} N et Ft =  {calcul_forces(0.0001)[1]} N."
            f"\nPour l'angle maximum immaginable de l'effort du sondeur sur le premier élement de cable "
            f"(la trainée s'opposerai à l'inclinaison d'où angle max) "
            f"on aurait \nFn = {calcul_forces(10)[0]} N et Ft = {calcul_forces(10)[1]} N"
            f"\nOn peud encadrer au moins par ces angles l'inclinaison du premier element de câble")

    #region --- PARTIE 1 : TRAÇAGE DÉFORMÉE ET LAYBACK UNIQUE ---

    # réglage les paramètres de la déformée
    Longueur_du_cable = 200
    Vitesse_du_bateau_en_noeuds = 3
    Nombre_de_section_de_la_deformee = 400

    if afficher_deformee:
        x, z = resoudre_deformee(Longueur_du_cable, Vitesse_du_bateau_en_noeuds, Nombre_de_section_de_la_deformee)
        tracer_deformee(x, z)

        layback = abs(x[0])
        profondeur = abs(z[0])
        print(f"\n--- RÉSULTATS ---")
        print(f"Layback (distance horizontale sonar-bateau) : {layback:.2f} m")
        print(f"Profondeur du sonar : {profondeur:.2f} m")
    #endregion

    #region --- PARTIE 2 : CALCUL ET TRAÇAGE COMPARATIF ---
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

    for i, L in enumerate(longueurs):
        for j, v in enumerate(vitesses):
            x, z = resoudre_deformee(L_tot=L, v_kt=v)
            laybacks[i, j] = abs(x[0])
            profondeurs[i, j] = abs(z[0])
    #endregion

    #region --- PARTIE 3 : AFFICHAGE DES COMPARAISONS ---

    if afficher_graphiques_comparaison:
        df_prof = pd.DataFrame(profondeur_constructeur, index=longueurs, columns=vitesses)
        df_layb = pd.DataFrame(layback_constructeur, index=longueurs, columns=vitesses)

        plt.figure(figsize=(10, 6))
        for vitesse in df_prof.columns:
            plt.plot(df_prof.index, df_prof[vitesse], marker='o', label=f"{vitesse} kt")
        plt.title("Profondeur (constructeur)")
        plt.xlabel("Longueur de câble (m)")
        plt.ylabel("Profondeur (m)")
        plt.grid(True)
        plt.legend(title="Vitesse")
        plt.tight_layout()
        plt.show()

        plt.figure(figsize=(10, 6))
        for vitesse in df_layb.columns:
            plt.plot(df_layb.index, df_layb[vitesse], marker='s', label=f"{vitesse} kt")
        plt.title("Layback (constructeur)")
        plt.xlabel("Longueur de câble (m)")
        plt.ylabel("Layback (m)")
        plt.grid(True)
        plt.legend(title="Vitesse")
        plt.tight_layout()
        plt.show()
    #endregion

    #region --- PARTIE 4 : TABLEAUX D'ÉCARTS RELATIFS ---

    if afficher_tableaux_ecarts:

        #region --- Tableau Layback ---
        df_laybacks = pd.DataFrame(laybacks, index=longueurs, columns=vitesses)
        print("\nLAYBACK (m) :")
        print(df_laybacks.round(3))

        df_laybacks_ref = pd.DataFrame(layback_constructeur, index=longueurs, columns=vitesses)
        print("\nLAYBACK CONSTRUCTEUR (m) :")
        print(df_laybacks_ref.round(3))

        ecarts_layback = laybacks - layback_constructeur
        df_ecart_layback = pd.DataFrame(ecarts_layback, index=longueurs, columns=vitesses)
        print("\nECARTS LAYBACK (m)) :")
        print(df_ecart_layback.round(3))

        laybacks_rel = (laybacks.T / longueurs)*100
        df_laybacks_rel = pd.DataFrame(laybacks_rel, index=longueurs, columns=vitesses)
        print("\nLAYBACK RELATIF % (Layback / Longueur câble) :")
        print(df_laybacks_rel.round(3))

        laybacks_ref_rel = (layback_constructeur.T / longueurs)*100
        df_laybacks_ref_rel = pd.DataFrame(laybacks_ref_rel, index=longueurs, columns=vitesses)
        print("\nLAYBACK RELATIF CONSTRUCTEUR % (Layback / Longueur câble) :")
        print(df_laybacks_ref_rel.round(3))

        ecarts_rel_layback = laybacks_rel - laybacks_ref_rel
        df_rel_layback = pd.DataFrame(ecarts_rel_layback, index=longueurs, columns=vitesses)
        print("\nECARTS LAYBACK RELATIFS % (Layback / Longueur câble) :")
        print(df_rel_layback.round(3))
        # endregion

        #region --- Tableau Prodondeur ---
        df_prof = pd.DataFrame(profondeurs, index=longueurs, columns=vitesses)
        print("\nPROFONDEUR (m) :")
        print(df_prof.round(3))

        df_prof_ref = pd.DataFrame(profondeur_constructeur, index=longueurs, columns=vitesses)
        print("\nPROFONDEUR CONSTRUCTEUR (m) :")
        print(df_prof_ref.round(3))

        ecarts_prof = profondeurs - profondeur_constructeur
        df_ecart_prof = pd.DataFrame(ecarts_prof, index=longueurs, columns=vitesses)
        print("\nECARTS PROFONDEUR (m)) :")
        print(df_ecart_prof.round(3))

        prof_rel = (profondeurs.T / longueurs)*100
        df_prof_rel = pd.DataFrame(prof_rel, index=longueurs, columns=vitesses)
        print("\nPROFONDEUR RELATIF % (Profondeur / Longueur câble) :")
        print(df_prof_rel.round(3))

        prof_ref_rel = (profondeur_constructeur.T / longueurs)*100
        df_prof_ref_rel = pd.DataFrame(prof_ref_rel, index=longueurs, columns=vitesses)
        print("\nPROFONDEUR RELATIF CONSTRUCTEUR % (Profondeur / Longueur câble) :")
        print(df_prof_ref_rel.round(3))

        ecarts_rel_prof = prof_rel - prof_ref_rel
        df_rel_prof = pd.DataFrame(ecarts_rel_prof, index=longueurs, columns=vitesses)
        print("\nECARTS PROFONDEUR RELATIFS % (Profondeur / Longueur câble) :")
        print(df_rel_prof.round(3))
        # endregion

        #region --- figures ---
        if afficher_graph:
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
        #endregion
    #endregion