import numpy as np
from simpson import integration


def main():
    A = 238
    Z = 92
    div = 10000000
    energy = 5.0 * (1.6 * 10**-13)
    # mass = 3727.3794118 * (1.6 * 10**-19)  # mass fo alpha particle in MeV
    mass = 6.6446e-27
    upper = 8 * np.power(A, (1 / 3)) * 10**-15
    lower = 1.07 * np.power(A, (1 / 3)) * 10**-15
    value = integration(transtition, upper, lower, div, Z, energy, mass)
    prob = np.exp(-2 * value)
    print(f"The transtition probabilit is:- {prob}")


def transtition(radius: float, Z: int, energy: float, mass: float):
    potential = col(radius, Z)
    temp = 2 * mass
    temp = temp / ((1.05e-34) ** 2)
    temp = temp * (potential - energy)
    temp = np.sqrt(temp)
    return temp


def col(r: float, Z: int) -> float:
    potential = (9 * 10**9) * 2 * Z * (1.6 * 10**-19) ** 2 / (r)
    # potential = potential / (1.6 * 10**-19)
    return potential


if __name__ == "__main__":
    main()
