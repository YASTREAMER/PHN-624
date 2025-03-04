import numpy as np
import pandas as pd

from hermite import hermite
from factorial import factorial
from simpson import integration


def main():
    np.sqrt(4)
    n, m = 5, 5
    upper = 8
    lower = -8
    KE = []
    total = []
    pot = []
    for i in range(n):
        temp = []
        for j in range(m):
            temp.append(integration(ke, upper, lower, 1000, j, i))

        KE.append(temp)
    for i in range(n):
        temp = []
        for j in range(m):
            temp.append(integration(potential, upper, lower, 1000, j, i))

        pot.append(temp)
    for i in range(n):
        temp = []
        for j in range(m):
            temp.append(KE[i][j] + pot[i][j])

        total.append(temp)

    kine = pd.DataFrame(KE)
    poten = pd.DataFrame(pot)
    totalEnergy = pd.DataFrame(total)
    print(f"The kinetic energy matrix is \n{ kine }")
    print(f"The potential energy matrix is \n{ poten}")
    print(f"The total energy matrix is \n{ totalEnergy}")


def ke(m: int, n: int, x: float) -> float:

    # if n - 2 < 0 or n - 1 < 0: return 0

    if m == n:
        return 0

    if n - 2 < 0:
        w1 = 0

    else:
        w1 = 2 * np.sqrt(n * n - 1) * wavefunction(n - 2, x)

    if n - 1 < 0:
        w2 = 0

    else:
        w2 = np.sqrt(8 * n) * x * wavefunction(n - 1, x)

    w3 = (x**2 - 1) * wavefunction(n, x)

    psiDP = w1 - w2 + w3

    energy = wavefunction(m, x) * psiDP * 0.5 * (1.0)
    return energy


def potential(m: int, n: int, x: float) -> float:
    energy = wavefunction(m, x) * ((x**2) / 2) * wavefunction(n, x)

    return energy


def wavefunction(n, x) -> float:

    if n < 0:
        return 0

    psi = np.sqrt(1 / (2**n * np.sqrt(np.pi) * factorial(n)))
    psi = psi * np.exp(-(x**2) / 2) * hermite(n, x)
    return psi


if __name__ == "__main__":
    main()
