import numpy as np


def main():

    # -------- Constants ---------

    Z1 = 84
    Z2 = 2
    A1 = 214
    A2 = 4
    NucleonMass = 938
    e = 1e-6
    hcut = 197.3
    Rnot = 1.2  # in fm
    energy = 8.8  # MeV
    energy = energy     

    vel = np.sqrt(2 * energy / A2 * NucleonMass)
    freq = Frequency(vel, Rnot, A1)
    rMass = ReducedMass(A1 * NucleonMass, A2 * NucleonMass)

    # -------- Calculation --------

    gFactor = GamowFactor(Z2, e, hcut, rMass, freq, vel=vel)
    print(f"The Gamow Factor is {gFactor}")

    lam = lambdaFunction(freq, gFactor)
    print(f"Lambda Value is {lam}")

    halfLife = 0.693 / lam
    print(f"The half life is {halfLife} seconds")


def GamowFactor(
    Zdaughter: int, e: float, hcut: float, rMass: float, freq: float, vel: float
) -> float:

    gFactor = 2 * (np.pi) * Zdaughter * e**2 / (hcut * vel)
    gFactor = gFactor - ((4 * e / hcut) * np.sqrt(Zdaughter * freq))

    return gFactor


def lambdaFunction(freq, gFactor):
    lam = freq * np.exp(-2 * gFactor)
    return lam

    pass


# ------- UTILS --------


def Frequency(vel: float, Rnot: float, Ap: float) -> float:
    R = Rnot * Ap ** (1 / 3)
    return vel / R


def ReducedMass(mass1: float, mass2: float) -> float:

    return mass1 * mass2 / (mass1 + mass2)


if __name__ == "__main__":
    main()
