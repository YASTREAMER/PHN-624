import numpy as np
import matplotlib.pyplot as plt


def main():
    energy = []
    mass = []

    for A in range(1, 251):

        if A % 2 == 0:
            Z = A / 2

        else:
            Z = (A - 1) / 2

        even = True if (A - Z) % 2 == 0 else False

        en = empherical(A, Z, even)

        energy.append((en / A))
        mass.append(A)
    print(energy[14])
    plt.plot(mass, energy)
    plt.show()


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


if __name__ == "__main__":
    main()
