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
D = 0.0063         # Diamètre du câble (m)
S = (D**2) * np.pi / 4 # section du cable
mass_lin = 0.5     # Masse linéique du câble (kg/m)
v_kt = 3          # Vitesse en noeud
v_bateau = v_kt * 0.5144  # Vitesse en m/s

# --- SONAR ---
m_sondeur_imm = 6.7             # Masse apparente (plongée)
D_sondeur = 0.06
S_sondeur = (D_sondeur**2) * np.pi / 4  # Surface frontale


# --- COEFFICIENTS HYDRODYNAMIQUES ---
Cd_sondeur = 0.7
Cd = 1.2     # Traînée normale
Cf = 0.01    # Traînée tangentielle

# --- 4. FORCES ----

#force sur l_element
def calcul_forces(theta):
    delta_s = L_tot / n
    v_n = v_bateau * np.cos(theta)
    v_t = v_bateau * np.sin(theta)

    F_n = - 0.5 * rho_eau * v_n**2 * Cd * D * delta_s
    F_t = - 0.5 * rho_eau * v_t**2 * Cf * S
    P = - (mass_lin - rho_eau * (np.pi * D**2 / 4)) * g * delta_s  # Poids apparent



    return F_n, F_t, P

def effort_sondeur():
    Fn = m_sondeur_imm * g
    Ft = 0.5 * rho_eau * Cd_sondeur * S_sondeur * v_bateau**2
    return Fn, Ft
print(effort_sondeur(), np+arctan(effort_sondeur()[0], effort_sondeur()[1]))



# --- 5. CALCUL DE LA DÉFORMÉE ---

def resoudre_deformee():
    theta = np.zeros(n + 1)
    x = np.zeros(n + 1)
    z = np.zeros(n + 1)
    eff_sond = (effort_sondeur())
    theta[0] = np.arctan(eff_sond[0], eff_sond[1])
    delta_s = L_tot / n

    F = [(eff_sond)]

    F_n, F_t, P = calcul_forces(theta[0])
    Fx = F_n * np.cos(theta[0]) - F_t * np.sin(theta[0])
    Fz = -P + F_n * np.sin(theta[0]) + F_t * np.cos(theta[0])
    F.append()
    F = (Fx, )

    for i in range(n):
        F_n, F_t, P = calcul_forces(theta[i])
        Fx = F_n * np.cos(theta[i]) - F_t * np.sin(theta[i])
        Fz = -P + F_n * np.sin(theta[i]) + F_t * np.cos(theta[i])



        theta[i+1] = np.arctan(Fz/Fx)
        x[i] = x[i] + delta_s * np.cos(theta[i+1])
        z[i] = z[i] + delta_s * np.sin(theta[i+1])

    return x, z


# --- 6. AFFICHAGE ---

def tracer_deformee(x, z):
    plt.figure(figsize=(10, 6))
    plt.plot(x, -z)  # z orienté vers le bas
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
