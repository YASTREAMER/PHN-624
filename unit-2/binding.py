import numpy as np
import matplotlib.pyplot as plt


def main():
    energy = []
    mass = []

    for A in range(2, 251):
        if A % 2 == 0:
            even = True
        else:
            even = False

        Z = stable(A)
        Z = int(round(Z))
        N = A - Z

        energy.append(empherical(A, Z, even)/N)
        mass.append(A)
    plt.plot(mass, energy)

    plt.show()


def empherical(A: int, Z: float, even: bool) -> float:

    volTerm = 15.75
    surfaceTerm = 17.8
    asymm = 23.7
    colTerm = 0.71
    pairTerm = 33.5

    if not (even):
        pairing = 0
    else:
        pairing = pairTerm / np.power(A, (3 / 4))

    if Z % 2 == 0:
        pairing = -pairing

    energy = (
        volTerm * A
        - surfaceTerm * np.power(A, (2 / 3))
        - asymm * np.power((A - Z * 2), 2) / A
        - colTerm * np.power(Z, 2) / np.power(A, (1 / 3))
        - pairing
    )

    return energy


def stable(A):
    asymm = 23.70
    colTerm = 0.71

    Z = (2 * asymm) / (2 * asymm + colTerm * np.power(A, (2 / 3)))
    Z = Z * A

    return Z


if __name__ == "__main__":
    main()
