import matplotlib.pyplot as plt
import numpy as np


def main():
    volTerm = 16
    surfaceTerm = 20
    asymm = 22.4
    colTerm = 0.751

    ZBeta = []  # Z for beta stablility line
    NBeta = []

    Zproton = []
    NProton = []

    NNeutron = []
    ZNeutron = []

    for A in range(1, 251):

        beta = betastable(A, colTerm, asymm)
        pro = protondrip(A, colTerm, asymm, volTerm, surfaceTerm)  # Value of Z
        neu = neutrondrip(A, volTerm, surfaceTerm, asymm, colTerm)  # Value of N

        # Beta stablility line
        Z = (A - beta) / 2
        N = A - Z
        # print(f"The value of N is {N} and the value of Z is {Z}")
        ZBeta.append(Z)
        NBeta.append(N)

        # Proton Drip line
        N = A - pro
        NProton.append(N)
        Zproton.append(pro)
        # print(f"The value of N is {N} and the value of Z is {pro}")

        # Neutron Drip line
        Z = A - neu
        NNeutron.append(neu)
        ZNeutron.append(Z)
        # print(f"The value of N is {Z} and the value of Z is {neu}")


    plt.plot(NBeta,ZBeta)
    plt.plot(NProton,Zproton)
    plt.plot(NNeutron,ZNeutron)
    plt.legend(["Beta stablility line ", "Proton drip line ", "Neutron drip line "])
    plt.xlabel("Z no of neutron")
    plt.ylabel("N no of proton")
    plt.show()


def coefalpha(A, volTerm, surfaceTerm, colTerm, asymm) -> float:

    alpha = (
        volTerm
        + 3 * asymm
        - (2 / 3) * surfaceTerm / (A ** (1 / 3))
        + colTerm * A ** (2 / 3) / 3
    )

    alpha = alpha / ((colTerm / 3 * A ** (2 / 3)) + 4 * asymm)

    return alpha


def coefbeta(A, colTerm, asymm) -> float:
    coef = 2 * colTerm * A ** (2 / 3) / (colTerm * A ** (2 / 3) + 12 * asymm)
    return coef


def coefgamma(A, volTerm, surfaceTerm, asymm, colTerm) -> float:
    coef = (volTerm + 3 * asymm - (2 / 3) * surfaceTerm * A ** (-1 / 3)) / (
        (1 / 3) * colTerm * A ** (2 / 3) + 4 * asymm
    )
    return coef


def betastable(A: int, colTerm: float, asymm: float, Mh=939.0, Mn=940.6) -> float:
    """

    Calculates the beta stablility line. We are taking mass of hydrogen neutron in MeV/c^2, thus we can take c = 1

    """

    temp = colTerm * np.power(A, (2 / 3)) - (Mn - Mh) * 1  # We are taking c = 1

    temp = temp / ((asymm * 4 / A) + (colTerm / np.power(A, (1 / 3))))

    return temp


def protondrip(A, colterm, asymm, volTerm, surfaceTerm) -> float:

    beta = coefbeta(A, colterm, asymm)
    gamma = coefgamma(A, volTerm, surfaceTerm, asymm, colterm)

    proton = A * ((1 + beta) - np.sqrt((1 + beta) ** 2 - gamma))
    return proton


def neutrondrip(A, volTerm, surfaceTerm, asymm, colTerm):

    alpha = coefalpha(A, volTerm, surfaceTerm, colTerm, asymm)

    neutron = A * (1 - np.sqrt(1 - alpha))
    return neutron


if __name__ == "__main__":
    main()
