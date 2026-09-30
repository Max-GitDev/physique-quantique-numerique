#!/usr/bin/env python3
"""
Résonance paramagnétique électronique (RPE) dans le rubis.

Diagonalisation de l'hamiltonien de spin S = 3/2 de l'ion Cr3+ (champ cristallin + effet Zeeman),
recherche des champs de résonance avec une cavité à 9,19 GHz, comparaison à des mesures
expérimentales, et calcul du couplage collectif spin-cavité Omega(T).

Auteur : Maxime Rimbaud - ESPCI Paris - PSL, TP numérique de physique quantique, option (2026).
"""

from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt
from scipy.constants import h, k as k_B, N_A

plt.rcParams.update({'font.size': 11, 'figure.dpi': 150})

FIG = Path(__file__).parent / "figures"
FIG.mkdir(exist_ok=True)


def sauvegarder(nom):
    """Sauvegarde la figure en PNG dans le dossier figures/."""
    plt.savefig(FIG / nom.replace('.pdf', '.png'), bbox_inches='tight', dpi=200)

# -- Parametres physiques -------------------------------------------------
D      = 5.73e9      # couplage champ cristallin (Hz)
g      = 1.98        # facteur de Lande Cr3+ dans le rubis
mu_B   = 9.274e-24   # magneton de Bohr (J/T)
f_cav  = 9.19e9      # frequence de resonance de la cavite (Hz)
b_y    = 100e-12     # amplitude du champ micro-onde (T)
rho    = 3.98e6      # masse volumique rubis (g/m^3)
M      = 102.0       # masse molaire Al2O3 (g/mol)
V      = 1e-9        # volume de l'echantillon (m^3)
N_ions = 0.001 * (2 * rho * V / M) * N_A   # nombre d'ions Cr3+

# -- Matrices de spin S=3/2 (sans dimension, en unites de hbar) -----------
Sz = np.diag([3/2, 1/2, -1/2, -3/2]).astype(complex)
Sp = np.zeros((4, 4), dtype=complex)
Sp[0, 1] = np.sqrt(3); Sp[1, 2] = 2.0; Sp[2, 3] = np.sqrt(3)
Sm = Sp.conj().T
Sx = (Sp + Sm) / 2.0
Sy = (Sp - Sm) / (2j)


def hamiltonien(B0_G, theta_deg):
    """H/h en Hz, B0 en Gauss, theta en degres."""
    B0       = B0_G * 1e-4
    theta    = np.deg2rad(theta_deg)
    H_crist  = -D * (Sz @ Sz - (5/4) * np.eye(4))   # terme champ cristallin
    nu_Z     = g * mu_B * B0 / h                      # frequence de Larmor
    H_zeeman = nu_Z * (np.cos(theta) * Sz + np.sin(theta) * Sx)
    return H_crist + H_zeeman


def diagonalise(B0_G, theta_deg):
    niveaux, etats = np.linalg.eigh(hamiltonien(B0_G, theta_deg))
    return niveaux.real, etats


# -- Etats propres aux points cles ----------------------------------------
print("=== Etats propres a B0=0, theta=40 deg ===")
niveaux0, etats0 = diagonalise(0.0, 40.0)
base = ['|3/2,+3/2>', '|3/2,+1/2>', '|3/2,-1/2>', '|3/2,-3/2>']
for k in range(4):
    coefficients = etats0[:, k]
    expression = " + ".join(f"({c.real:+.3f}){b}"
                            for c, b in zip(coefficients, base) if abs(c) > 0.01)
    print(f"  E{k}/h = {niveaux0[k]/1e9:.3f} GHz :  {expression}")

print("\n=== Etats propres a B0=6000G, theta=40 deg ===")
niveaux6, etats6 = diagonalise(6000.0, 40.0)
for k in range(4):
    coefficients = etats6[:, k]
    termes = [(abs(c), c, b) for c, b in zip(coefficients, base)]
    expression = " + ".join(f"({c.real:+.3f}){b}"
                            for _, c, b in sorted(termes, reverse=True)[:4]
                            if abs(c) > 0.01)
    print(f"  E{k}/h = {niveaux6[k]/1e9:.3f} GHz :  {expression}")


# -- Recherche des champs de resonance : |E_i - E_j| = f_cav -------------
def cherche_resonances(theta_deg, B_min=0, B_max=6000, n_points=3000):
    """Retourne la liste (i, j, B_res) des paires en resonance avec la cavite."""
    B_tab     = np.linspace(B_min, B_max, n_points)
    E_tab     = np.array([diagonalise(B, theta_deg)[0] for B in B_tab])
    resultats = []
    for i in range(4):
        for j in range(i + 1, 4):
            ecart = np.abs(E_tab[:, i] - E_tab[:, j]) - f_cav
            for k in np.where(np.diff(np.sign(ecart)))[0]:
                if ecart[k + 1] != ecart[k]:
                    B_res = B_tab[k] - ecart[k] * (B_tab[k+1] - B_tab[k]) / (ecart[k+1] - ecart[k])
                    resultats.append((i, j, float(B_res)))
    return resultats


# -- Figure 1 : niveaux d'energie a theta=40 deg --------------------------
theta_ref = 40.0
B_tab1    = np.linspace(0, 6000, 600)
E_tab1    = np.array([diagonalise(B, theta_ref)[0] for B in B_tab1])
res_40    = cherche_resonances(theta_ref)

print(f"\n=== Transitions a theta=40 deg, f_cav={f_cav/1e9} GHz ===")
for (i, j, B_res) in res_40:
    print(f"  niveaux {i+1}-{j+1} : B0 = {B_res:.0f} G")

fig, ax = plt.subplots(figsize=(7, 5))
couleurs = ['#1d7874', '#ee2e31', '#9b59b6', '#f0a500']
for k in range(4):
    ax.plot(B_tab1, E_tab1[:, k] / 1e9, color=couleurs[k], lw=1.5, label=f'niveau {k+1}')
ax.axhline(f_cav / 1e9, color='k', ls='--', lw=1,
           label=f'$f_{{\\rm cav}} = {f_cav/1e9:.2f}$ GHz')
for (i, j, B_res) in res_40:
    ax.axvline(B_res, color='#888888', ls=':', lw=0.8)
ax.set_xlabel('$B_0$ (G)')
ax.set_ylabel('$E/h$ (GHz)')
ax.set_title(f'Niveaux d\'energie, $\\theta = {theta_ref:.0f}^{{\\circ}}$')
ax.legend(fontsize=9, loc='upper left')
plt.tight_layout()
sauvegarder('fig1_niveaux.pdf')
plt.close()

# -- Figure 1b : frequences de transition a theta=40 deg ------------------
paires          = [(0,1),(0,2),(0,3),(1,2),(1,3),(2,3)]
etiquettes      = ['niveaux 1-2', 'niveaux 1-3', 'niveaux 1-4',
                   'niveaux 2-3', 'niveaux 2-4', 'niveaux 3-4']
couleurs_paires = ['#1d7874', '#ee2e31', '#9b59b6', '#f0a500', '#2c7bb6', '#a0522d']

fig, ax = plt.subplots(figsize=(7, 5))
for (i, j), etiquette, couleur in zip(paires, etiquettes, couleurs_paires):
    freq_diff = np.abs(E_tab1[:, i] - E_tab1[:, j]) / 1e9
    ax.plot(B_tab1, freq_diff, color=couleur, lw=1.5, label=etiquette)
ax.axhline(f_cav / 1e9, color='k', ls='--', lw=1.2,
           label=f'$f_{{\\rm cav}} = {f_cav/1e9:.2f}$ GHz')
for (i, j, B_res) in res_40:
    ax.axvline(B_res, color='#888888', ls=':', lw=0.8)
ax.set_xlabel('$B_0$ (G)')
ax.set_ylabel('$|E_i - E_j|/h$ (GHz)')
ax.set_title(f'Frequences de transition, $\\theta = {theta_ref:.0f}^{{\\circ}}$')
ax.set_ylim(0, 30)
ax.legend(fontsize=9, loc='upper right')
plt.tight_layout()
sauvegarder('fig1b_freq.pdf')
plt.close()

# -- Figures 2 & 3 : carte des resonances dans le plan (theta, B0) --------
angles_tab = np.linspace(0, 190, 192)
angles_res, champs_res = [], []
for theta in angles_tab:
    for (i, j, B_res) in cherche_resonances(theta):
        angles_res.append(theta)
        champs_res.append(B_res)

donnees_exp = {   # mesures de la promotion 2018 fournies avec l'énoncé
    0:  [750, 1050, 3700, 3300],
    20: [1050, 2620, 4220, 5750],
    40: [460, 1500, 2600, 4400],
    60: [550, 2040, 3600],
    80: [900, 2140, 4820],
}

fig, ax = plt.subplots(figsize=(8, 5))
ax.scatter(angles_res, champs_res, s=2, color='#2c7bb6', alpha=0.6, label='Simulation')
ax.set_xlabel('Angle $\\theta$ ($^{\\circ}$)')
ax.set_ylabel('$B_0$ (G)')
ax.set_title('Transitions RPE : $|E_i - E_j| = hf_{\\rm cav}$')
ax.set_xlim(0, 190); ax.set_ylim(0, 6000)
ax.legend(); plt.tight_layout()
sauvegarder('fig2_transitions.pdf')
plt.close()

fig, ax = plt.subplots(figsize=(8, 5))
ax.scatter(angles_res, champs_res, s=2, color='#2c7bb6', alpha=0.5, label='Simulation')
premier = True
for theta_exp, champs_exp in donnees_exp.items():
    ax.scatter([theta_exp] * len(champs_exp), champs_exp, color='#d7191c', s=60,
               zorder=5, marker='D', linewidths=1.2,
               label='Donnees experimentales' if premier else '')
    premier = False
ax.set_xlabel('Angle $\\theta$ ($^{\\circ}$)')
ax.set_ylabel('$B_0$ (G)')
ax.set_title('Simulation vs donnees experimentales')
ax.set_xlim(0, 190); ax.set_ylim(0, 6000)
ax.legend(); plt.tight_layout()
sauvegarder('fig3_exp.pdf')
plt.close()


# -- Couplage spin-cavite Omega(T) ----------------------------------------
temperatures = np.logspace(np.log10(0.01), np.log10(30), 400)   # 10 mK -> 30 K


def couplage_omega(theta_deg, B0_G, i, j, temperatures):
    """Retourne Omega(T)/h en Hz pour la transition (i, j) a (theta, B0)."""
    niveaux, etats = diagonalise(B0_G, theta_deg)
    mel       = abs(etats[:, i].conj() @ (Sy @ etats[:, j]))  # <phi_i|Sy|phi_j>
    E_decalee = (niveaux - niveaux[0]) * h                    # energies en J

    def diff_populations(T):
        poids  = np.exp(-E_decalee / (k_B * T))
        Z_part = poids.sum()                                   # fonction de partition
        return (poids[i] - poids[j]) / Z_part

    polarisation = np.array([diff_populations(T) for T in temperatures])
    return (g * mu_B / h) * b_y * np.sqrt(N_ions * np.abs(polarisation)) * mel


PALETTE_COUPLAGE = ['#e63946', '#457b9d', '#2a9d8f', '#e9c46a']

def trace_couplage(theta_deg, nom_fichier, titre):
    transitions = cherche_resonances(theta_deg)
    print(f"\n=== Omega(T) a theta={theta_deg} deg ===")
    fig, ax = plt.subplots(figsize=(7, 5))
    for k, (i, j, B_res) in enumerate(transitions):
        omega     = couplage_omega(theta_deg, B_res, i, j, temperatures)
        etiquette = f'niveaux {i+1}-{j+1}, $B_0={B_res:.0f}$ G'
        ax.plot(temperatures, omega / 1e6, lw=1.5,
                color=PALETTE_COUPLAGE[k % len(PALETTE_COUPLAGE)], label=etiquette)
        print(f"  {etiquette}  |  Omega_max/h = {omega.max()/1e6:.2f} MHz")
    ax.set_xscale('log')
    ax.set_xlabel('$T$ (K)')
    ax.set_ylabel('$\\Omega / h$ (MHz)')
    ax.set_title(titre)
    ax.legend(fontsize=8)
    plt.tight_layout()
    sauvegarder(nom_fichier)
    plt.close()


trace_couplage(36.0, 'fig4_omega36.pdf', 'Couplage spin-cavite, $\\theta = 36^{\\circ}$')
trace_couplage(98.0, 'fig5_omega98.pdf', 'Couplage spin-cavite, $\\theta = 98^{\\circ}$')

print(f"\nN_ions = {N_ions:.3e}")
print("Figures sauvegardees.")
