import numpy as np
import os 


def main():
    """

    This function calculates the binding energy per nucleon

    """

    volTerm = 15.75
    surfaceTerm = 17.8
    asymm = 23.7
    colTerm = 0.71
    A = [14, 26, 25, 110, 202]
    Z = [7, 12, 17, 48, 80]
    Nuclei = ["N", "Mg", "Cl", "Cd", "Hg"]

    if os.path.isfile("Question2.txt"):
        os.remove("Question2.txt")


    file = open("Question2.txt", "a")

    for i in range(A.__len__()):
        energy = calculate(A[i], Z[i], Nuclei[i])
        val = protondrip(A[i], colTerm, asymm, volTerm, surfaceTerm)
        file.write(f"The proton drip line of {A[i]} is {val}\n")
        file.write(f"The binding energy per nucleon of {Nuclei[i]} is : {energy/A[i]} \n")


        print(f"The proton drip line of {A[i]} is {val}\n")
        print(f"The binding energy per nucleon of {Nuclei[i]} is : {energy/A[i]} \n")


def calculate(A: int, Z: int, Nuclei: str) -> float:
    even = True if (A - Z) % 2 == 0 else False
    en = empherical(A, Z, even)
    return en / A


def empherical(A: int, Z: float, even: bool) -> float:

    volTerm = 15.75
    surfaceTerm = 17.8
    asymm = 23.7
    colTerm = 0.71
    pairTerm = 33.5

    if even:
        pairing = pairTerm / np.power(A, (3 / 4))
    else:
        pairing = 0

    if Z % 2 != 0 and even:
        pairing = -pairing

    sur = surfaceTerm * np.power(A, (2 / 3))
    vol = volTerm * A

    aterm = np.power((A - 2 * Z), 2) / A
    aterm = aterm * asymm

    col = colTerm * np.power(Z, 2) / np.power(A, (1 / 3))

    energy = vol - sur - aterm - col + pairing
    return energy


def protondrip(A, colterm, asymm, volTerm, surfaceTerm) -> tuple:

    beta = coefbeta(A, colterm, asymm)
    gamma = coefgamma(A, volTerm, surfaceTerm, asymm, colterm)

    proton1 = A * ((1 + beta) - np.sqrt((1 + beta) ** 2 - gamma))
    proton2 = A * ((1 + beta) + np.sqrt((1 + beta) ** 2 - gamma))
    return proton1, proton2


def coefbeta(A, colTerm, asymm) -> float:
    coef = 2 * colTerm * A ** (2 / 3) / (colTerm * A ** (2 / 3) + 12 * asymm)
    return coef


def coefgamma(A, volTerm, surfaceTerm, asymm, colTerm) -> float:
    coef = (volTerm + 3 * asymm - (2 / 3) * surfaceTerm * A ** (-1 / 3)) / (
        (1 / 3) * colTerm * A ** (2 / 3) + 4 * asymm
    )
    return coef


if __name__ == "__main__":
    main()
