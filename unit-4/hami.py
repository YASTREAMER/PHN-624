import numpy as np
import matplotlib.pyplot as plt


def main(num: int = 100) -> None:
    N = int(input("Enter the ending value of N \n"))
    deltaosc = np.linspace(-1, 1, num=num)
    energy = []
    for i in range(1, N + 1):
        for nz in range(1, i + 1):
            energy = Energyhw(i, nz, deltaosc)
            plt.plot(deltaosc, energy)

    plt.title("Anisotropic Oscillator with hw = 1")
    plt.show()

    for i in range(1, N + 1):
        for nz in range(1, i + 1):
            energy = Energy(i, nz, deltaosc)
            plt.plot(deltaosc, energy)

    plt.title("Anisotropic Oscillator with hw != 1")
    plt.show()


def Energy(N: int, nz: int, deltaosc: np.ndarray) -> list:
    energy = []
    for delta in deltaosc:
        hw = HcutOmega(delta)
        energy.append(((N + 1.5) - (1 / 3 * delta * (3 * nz - N))) * hw)
    return energy


def Energyhw(N: int, nz: int, deltaosc: np.ndarray) -> list:
    energy = []
    hw = 1
    for delta in deltaosc:
        energy.append(((N + 1.5) - (1 / 3 * delta * (3 * nz - N))) * hw)
    return energy


def HcutOmega(delta: float, A: int = 80) -> float:
    deltafunc = (1 + 2 / 3 * delta) ** 2
    deltafunc = deltafunc * (1 - (4 / 3) * delta)
    deltafunc = deltafunc ** (-1 / 6)
    hw = 41 * A ** (-1 / 3) * deltafunc
    return hw


if __name__ == "__main__":
    main(num=1000)
