import numpy as np
import matplotlib.pyplot as plt


def main():

    for A in range(1, 251, 1):
        if A == 1:
            Z = 1
            even = False
        elif A % 2 == 0:
            Z = A / 2
            even = True
        else:
            Z = (A - 1) / 2
            even = False

        energy = empherical(A, Z, even)
        plt.plot((energy/A, A))

    plt.show()


def empherical(A: int, Z: float, even: bool) -> float:

    volTerm = 157.716
    surfaceTerm = 414.976
    asymm = 4612.42
    colTerm = 14.131
    pairTerm = 13293.032

    if not (even):
        pairing = 0
    else:
        pairing = pairTerm / np.power(A, (3 / 4))

    if Z % 2 == 0:
        pairing = -pairing

    energy = (
        volTerm * A
        - surfaceTerm * np.power(A, (2 / 3))
        + asymm * np.power((A - 2 * Z), 2) / A
        - colTerm * np.power(Z, 2) / np.power(A, (1 / 3))
        - pairing
    )

    return energy


if __name__ == "__main__":
    main()
