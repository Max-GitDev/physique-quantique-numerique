#!/usr/bin/env python3
"""
États liés d'un puits quantique fini (hétérostructure AlGaAs / GaAs / AlGaAs).

Puits de largeur 2a et de profondeur V0, tronqué par un puits infini de largeur 2b.
On diagonalise l'hamiltonien dans la base des états propres du puits infini
phi_n(z) = sin(n pi (z + b) / 2b) / sqrt(b), puis on vérifie les conditions de
quantification |cos(q a)| = q / q0 (états pairs) et |sin(q a)| = q / q0 (états impairs).

Auteur : Maxime Rimbaud - ESPCI Paris - PSL, TP numérique de physique quantique (2025).
"""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from scipy.constants import electron_volt as eV, hbar, m_e
from scipy.integrate import quad

FIG = Path(__file__).parent / "figures"
FIG.mkdir(exist_ok=True)
plt.rcParams.update({"font.size": 11, "figure.dpi": 120})

a = 1e-9        # demi-largeur du puits fini (m)
b = 5 * a       # demi-largeur du puits infini (m)
V0 = 0.5        # profondeur du puits (eV)
N = 150         # taille de la base


def phi(n, z):
    return np.sin(n * np.pi * (z + b) / (2 * b)) / np.sqrt(b)


def energie_puits_infini(n):
    """Énergie du n-ième état du puits infini de largeur 2b, décalée de -V0 (eV)."""
    return (n * np.pi * hbar) ** 2 / (8 * m_e * b**2) / eV - V0


def hamiltonien(N):
    """H_nm = E_n delta_nm + int phi_n V phi_m, avec V = +V0 pour a < |z| < b (barrières)."""
    H = np.diag([energie_puits_infini(n) for n in range(1, N + 1)])
    for n in range(1, N + 1):
        for m in range(n, N + 1):
            if (n + m) % 2:                     # couplage nul entre états de parités différentes
                continue
            barrieres = sum(quad(lambda z: phi(n, z) * phi(m, z), zmin, zmax, limit=200)[0]
                            for zmin, zmax in ((-b, -a), (a, b)))
            H[n - 1, m - 1] += V0 * barrieres
            if m != n:
                H[m - 1, n - 1] = H[n - 1, m - 1]
    return H


if __name__ == "__main__":
    energies, vecteurs = np.linalg.eigh(hamiltonien(N))
    lies = np.where(energies < 0)[0]
    print(f"5 énergies propres les plus basses (eV) : {np.round(energies[:5], 4)}")
    print(f"Nombre d'états liés (E < 0) : {len(lies)}")

    z = np.linspace(-b, b, 1000)
    fonctions = []
    for k in lies:
        psi = sum(vecteurs[n - 1, k] * phi(n, z) for n in range(1, N + 1))
        psi /= np.sqrt(np.trapezoid(psi**2, z))
        fonctions.append(psi)

    # Vérification des conditions de quantification du puits fini (TD n°1)
    q0 = np.sqrt(2 * m_e * V0 * eV) / hbar
    print(f"\nq0 = {q0:.3e} m^-1")
    for k, E in zip(lies, energies[lies]):
        q = np.sqrt(2 * m_e * (E + V0) * eV) / hbar      # vecteur d'onde dans le puits
        pair = k % 2 == 0
        membre_gauche = abs(np.cos(q * a)) if pair else abs(np.sin(q * a))
        print(f"État {k} ({'pair' if pair else 'impair'}) : E = {E:+.4f} eV, "
              f"|{'cos' if pair else 'sin'}(qa)| = {membre_gauche:.3f}, q/q0 = {q / q0:.3f}")

    # Figure : potentiel, niveaux et fonctions d'onde (amplitude arbitraire)
    potentiel = np.where(np.abs(z) < a, -V0, 0.0)
    plt.figure(figsize=(7, 5))
    plt.plot(z * 1e9, potentiel, "k", lw=2, label="V(z)")
    echelle = 0.08 / max(np.abs(f).max() for f in fonctions)
    for k, psi in zip(lies, fonctions):
        plt.axhline(energies[k], color="grey", ls=":", lw=0.8)
        plt.plot(z * 1e9, energies[k] + echelle * psi, lw=2, label=f"état {k} : E = {energies[k]:.3f} eV")
    plt.xlim(-b * 1e9, b * 1e9)
    plt.ylim(-0.7, 0.2)
    plt.xlabel("z (nm)")
    plt.ylabel("Énergie (eV)")
    plt.title("États liés du puits fini GaAs (a = 1 nm, V0 = 0,5 eV)")
    plt.legend(loc="lower right", fontsize=9)
    plt.savefig(FIG / "puits_fini_etats_lies.png", bbox_inches="tight")
    plt.close()
    print(f"\nFigure enregistrée dans {FIG}")
