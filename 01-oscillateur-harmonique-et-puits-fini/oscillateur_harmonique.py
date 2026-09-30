#!/usr/bin/env python3
"""
Oscillateur harmonique quantique tronqué par un puits infini (méthode de Marsiglio).

On développe la fonction d'onde sur la base des états propres du puits infini
phi_n(x) = sqrt(2) sin(n pi x), 0 < x < 1, et on diagonalise la matrice de
l'hamiltonien normalisé  h_nm = n^2 delta_nm + 2 ∫ sin(n pi x) nu_OH(x) sin(m pi x) dx
avec nu_OH(x) = (pi^2 / 4) R^2 (x - 1/2)^2  et  R = hbar omega / E_1.

Référence : F. Marsiglio, "The harmonic oscillator in quantum mechanics: A third way",
Am. J. Phys. 77, 253 (2009).

Auteur : Maxime Rimbaud — ESPCI Paris - PSL, TP numérique de physique quantique (2025).
"""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import quad

FIG = Path(__file__).parent / "figures"
FIG.mkdir(exist_ok=True)
plt.rcParams.update({"font.size": 11, "figure.dpi": 120})

R = 30          # rapport hbar*omega / E_1 imposé par l'énoncé
N = 60          # taille de la base pour l'étude des énergies


def nu_oh(x, R):
    """Potentiel harmonique normalisé par E_1, centré en x = 1/2."""
    return (np.pi**2 / 4) * R**2 * (x - 0.5) ** 2


def hamiltonien(R, N):
    """Matrice h (N x N) calculée par intégration numérique."""
    h = np.zeros((N, N))
    for n in range(1, N + 1):
        for m in range(n, N + 1):
            integrale = quad(lambda x: np.sin(n * np.pi * x) * np.sin(m * np.pi * x) * nu_oh(x, R), 0, 1)[0]
            h[n - 1, m - 1] = h[m - 1, n - 1] = (n**2 if n == m else 0) + 2 * integrale
    return h


def hamiltonien_analytique(R, N):
    """Même matrice à partir des expressions analytiques démontrées dans le rapport."""
    n = np.arange(1, N + 1)
    nn, mm = np.meshgrid(n, n, indexing="ij")
    meme_parite = (nn - mm) % 2 == 0
    with np.errstate(divide="ignore", invalid="ignore"):
        h = np.where(meme_parite, 2 * R**2 * nn * mm / (nn**2 - mm**2) ** 2, 0.0)
    np.fill_diagonal(h, n**2 + R**2 / 8 * (np.pi**2 / 6 - 1 / n**2))
    return h


def etats_propres(R, N):
    """Énergies propres triées et vecteurs propres associés (en colonnes)."""
    energies, vecteurs = np.linalg.eigh(hamiltonien(R, N))
    return energies, vecteurs


def fonction_onde(coefficients, x):
    """Psi(x) = sum_j c_j sqrt(2) sin(j pi x)."""
    j = np.arange(1, len(coefficients) + 1)
    return np.sqrt(2) * np.einsum("ij,j->i", np.sin(np.outer(x, j) * np.pi), coefficients)


# --- Solutions analytiques de l'oscillateur harmonique (a = 1) -------------------------
def psi0_th(x, R):
    return (np.pi * R / 2) ** 0.25 * np.exp(-np.pi**2 * R * (x - 0.5) ** 2 / 4)


def psi1_th(x, R):
    return (np.pi**5 * R**3 / 2) ** 0.25 * (x - 0.5) * np.exp(-np.pi**2 * R * (x - 0.5) ** 2 / 4)


def psi2_th(x, R):
    return (np.pi * R / 8) ** 0.25 * (np.pi**2 * R * (x - 0.5) ** 2 - 1) * np.exp(-np.pi**2 * R * (x - 0.5) ** 2 / 4)


if __name__ == "__main__":
    # 1) Potentiel normalisé pour plusieurs valeurs de R
    x = np.linspace(0, 1, 1000)
    plt.figure()
    for r in (3, 6, 9, 21):
        plt.plot(x, nu_oh(x, r), label=f"R = {r}")
    plt.ylim(0, 10)
    plt.xlabel("x / a")
    plt.ylabel(r"$\nu_{OH}(x)$")
    plt.title("Potentiel harmonique normalisé dans le puits infini")
    plt.legend()
    plt.savefig(FIG / "potentiel_harmonique.png", bbox_inches="tight")
    plt.close()

    # 2) Matrice H pour N = 6 et vérification par les formules analytiques
    h6 = hamiltonien(R, 6)
    print("Matrice h normalisée (N = 6, R = 30) :")
    print(np.array2string(h6, precision=3, suppress_small=True))
    ecart = np.abs(hamiltonien(R, N) - hamiltonien_analytique(R, N)).max()
    print(f"Écart max. intégration numérique / formules analytiques (N = {N}) : {ecart:.1e}")

    # 3) Énergies propres normalisées et comparaison avec l'oscillateur harmonique
    eps, _ = etats_propres(R, N)
    n = np.arange(1, N + 1)
    eps_oh = R * (n - 1 + 0.5)                      # E_OH(k) / E_1 = R (k + 1/2), k = n - 1

    n0 = int(np.argmax(eps > nu_oh(1.0, R))) + 1    # premier niveau au-dessus du bord du puits
    regime_puits = n > n0
    C = np.mean(eps[regime_puits] - n[regime_puits] ** 2)   # moindres carrés pour f(n) = n^2 + C

    print(f"\n5 premières énergies normalisées : {np.round(eps[:5], 4)}")
    print(f"Oscillateur harmonique R(k+1/2)  : {eps_oh[:5]}")
    print(f"Transition de régime : nu_OH(1) = {nu_oh(1.0, R):.1f}, n0 = {n0}, eps(n0) = {eps[n0 - 1]:.1f}")
    print(f"Ajustement f(n) = n^2 + C sur n > n0 : C = {C:.1f}  (estimation diagonale pi^2 R^2 / 48 = {np.pi**2 * R**2 / 48:.1f})")

    plt.figure()
    plt.plot(n, eps, "o", ms=4, label=r"$\varepsilon(n)$ numérique")
    plt.plot(n, eps_oh, "-", label=r"oscillateur harmonique $R(k+\frac{1}{2})$")
    plt.plot(n, n**2 + C, "--", label=rf"$n^2 + C$, $C = {C:.0f}$")
    plt.axhline(nu_oh(1.0, R), color="grey", ls=":", label=r"bord du puits $\nu_{OH}(1)$")
    plt.axvline(n0, color="grey", ls=":")
    plt.xlabel("n")
    plt.ylabel(r"$\varepsilon(n) = E_n / E_1$")
    plt.title(f"Énergies propres normalisées (R = {R}, N = {N})")
    plt.legend()
    plt.savefig(FIG / "energies_propres.png", bbox_inches="tight")
    plt.close()

    # 4) Convergence des fonctions propres avec la taille de la base
    x = np.linspace(0, 1, 400)
    theoriques = (psi0_th, psi1_th, psi2_th)
    fig, axes = plt.subplots(1, 3, figsize=(14, 4), sharey=True)
    for etat, (ax, psi_th) in enumerate(zip(axes, theoriques)):
        reference = psi_th(x, R)
        for taille in (3, 5, 8, 15):
            _, vecteurs = etats_propres(R, taille)
            psi = fonction_onde(vecteurs[:, etat], x)
            psi *= np.sign(psi @ reference)          # phase globale arbitraire
            ax.plot(x, psi, label=f"N = {taille}")
        ax.plot(x, reference, "k--", lw=2, label="théorie")
        ax.set_title(rf"$\Psi_{etat}(x)$")
        ax.set_xlabel("x / a")
    axes[0].legend()
    fig.suptitle(f"Convergence des fonctions propres avec la taille de la base (R = {R})")
    fig.savefig(FIG / "fonctions_propres.png", bbox_inches="tight")
    plt.close(fig)

    print(f"\nFigures enregistrées dans {FIG}")
