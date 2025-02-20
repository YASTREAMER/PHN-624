import matplotlib.pyplot as plt
import numpy as np


def main():
    volTerm = 16
    surfaceTerm = 20
    asymm = 22.4
    colTerm = 0.751

    stable = []
    Zval = []
    Nval = []
    proton = []

    for A in range(1, 251):

        if A % 2 == 0:
            Z = A / 2

        else:
            Z = (A - 1) / 2

        Z = int(Z)

        N = A - Z

        beta = betastable(A, colTerm, asymm)
        pro = protondrip(A, Z, colTerm, asymm, volTerm, surfaceTerm)
        proton.append(pro)
        stable.append(beta)
        Zval.append(Z)
        Nval.append(N)

    plt.plot(Nval, stable)
    plt.plot(Nval, proton)
    plt.legend(["Beta stablility line ", "Proton stablility line "])
    plt.xlabel("N no of neutron")
    plt.ylabel("Z no of proton")
    plt.show()

def coefbeta(A, colterm, asymm) -> float:
    coef = 2 * colterm * A ** (2 / 3) / (colterm * A ** (2 / 3) + 12 * asymm)
    return coef


def coefgamma(A, volTerm, surfaceTerm, asymm, colterm) -> float:
    coef = (volTerm + 3 * asymm - (2 / 3) * surfaceTerm * A ** (-1 / 3)) / (
        (1 / 3) * colterm * A ** (2 / 3) + 4 * asymm
    )
    return coef


def betastable(A: int, colTerm: float, asymm: float, Mh=939.0, Mn=940.6) -> float:
    """

    Calculates the beta stablility line. We are taking mass of hydrogen neutron in MeV/c^2, thus we can take c = 1

    """

    temp = colTerm * np.power(A, (2 / 3)) - (Mn - Mh) * 1  # We are taking c = 1

    temp = temp / ((asymm * 4 / A) + (colTerm / np.power(A, (1 / 3))))

    return temp

def protondrip(A, Z, colterm, asymm, volTerm, surfaceTerm) -> float:

    beta = coefbeta(A, colterm, asymm)
    gamma = coefgamma(A, volTerm, surfaceTerm, asymm, colterm)

    proton = A * ((1 + beta) + np.sqrt((1 + beta) ** 2 + gamma))
    return proton

if __name__ == "__main__":
    main()
