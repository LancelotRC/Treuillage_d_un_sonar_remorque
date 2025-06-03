# Treuillage d'un Sonar Remorqué

Ce dépôt contient un projet académique réalisé en 2025, portant sur l'étude du treuillage d'un sonar latéral remorqué.  
Il inclut des analyses mécaniques, des modélisations numériques et des documents de synthèse.

---

## Contenu du dépôt

- `Déformée du cable.py` : Script Python pour la modélisation de la déformée du câble sous contraintes mécaniques.
- `auxil/` : Fonctions auxiliaires utilisées dans les scripts principaux.
- `LaTeX/` : Fichiers sources LaTeX pour la génération des documents PDF.
- `out/` : Dossier de sortie contenant les résultats des simulations et les figures générées.
- `Fiche_synthèse_sujet_Treuillage_sonar_remorqué_2025.pdf` : Document de synthèse du sujet traité.
- `Modélisation des mouillages.pdf` : Analyse des mouillages associés au système de treuillage.
- `UE2.2_Présentation_Sujet_Treuillage_sonar_remorqué_2025.pdf` : Présentation dans le cadre de l’UE2.2.
- `UE22_Treuillage_sonar_lateral_2025_vf.pdf` : Rapport final du projet.
- `sidescan.pdf` : Documentation technique sur le sonar latéral utilisé.

---

## Technologies utilisées

- **Python** : Pour la modélisation et les simulations.
- **LaTeX** : Pour la rédaction des rapports et présentations.

---

## Prérequis

- Python 3.x avec les bibliothèques suivantes :
  - `numpy`
  - `matplotlib`
  - `scipy` (si utilisé dans les scripts)
- Distribution LaTeX (TeX Live, MiKTeX…) pour compiler les documents PDF.

---

## Installation et utilisation

1. Cloner le dépôt :
   ```bash
   git clone https://github.com/LancelotRC/Treuillage_d_un_sonar_remorque.git

2. Installer les dépendances Python :

- A installer manuellement
    ```bash
    pip install numpy matplotlib

3. Exécuter les scripts Python pour générer les courbes ou visualiser les résultats :

python "Déformée du cable.py"

4. Compiler les documents LaTeX pour obtenir les rapports au format PDF (optionnel) :
    ```bash
    cd LaTeX
    pdflatex nom_du_document.tex
   
   ```

> 💡 Assurez-vous d’avoir une distribution LaTeX installée (TeX Live, MiKTeX…) pour compiler les documents.

---

## Auteurs

- Elouan S. — developpement informatique des pièces
- Thomas F. — Analyse et exploitation des données du sonar
- **pinsarew** — Calcul des expressions, modélisation et rédaction des documents.
- **LancelotRC** — Développement des scripts et rédaction des documents.


---

## Remarques

Ce projet a été réalisé dans un cadre pédagogique (année 2025).
L'utilisation ou la copie d'information issue de ce projet n'est pas permise pour des fins commerciaux.
Il serait apprécié de ne pas copier les informations présentes pour compléter son propre projet.