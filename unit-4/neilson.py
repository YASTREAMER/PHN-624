import numpy as np
from hami import *
from CGcoeff import cal


def main(sigma=[0.5, -0.5], num=100):
    delta = np.linspace(-0.3, 0.3, num)
    nuVals = [0.0, 0.0, 0.0, 0.35, 0.625, 0.63, 0.448, 0.434]
    k = 0.05
    value = []
    for N, val in enumerate(nuVals):
        for NPrime, valPrime in enumerate(nuVals):
            CalNillson(int(N), int(NPrime), sigma, sigma)


def CalNillson(N, NPrime, sigma, sigmaPrime):
    l = CalLval(N)
    lPrime = CalLval(NPrime)
    omegaVals = CalOmega(N)
    print(omegaVals)
    for sigmaVal in sigma:
        for sigmaVal in sigma:


def HCutDelta(
    delta: float, hcutomega0: float, l, lprime, lambdaVal, lambdaPrime
) -> float:
    energy = (
        delta
        * hcutomega0
        * (4 / 3)
        * np.sqrt(np.pi / 5)
        * SphericalHarmonics(l, lprime, lambdaVal, lambdaPrime)
    )
    return energy


def CalL2(N, l, lambdaVal, sigma, NPrime, lPrime, lambdaPrime, sigmaPrime) -> float:
    val = l * (l + 1)
    val = val * DiracDelta(N, NPrime)
    val = val * DiracDelta(l, lPrime)
    val = val * DiracDelta(lambdaVal, lambdaPrime)
    val = val * DiracDelta(sigma, sigmaPrime)
    return val


def CalLdotS(N, l, lambdaVal, sigma, NPrime, lPrime, lambdaPrime, sigmaPrime) -> float:

    val = 0.5 * np.sqrt((l - lambdaVal) * (l + lambdaPrime + 1))
    val1 = (
        val * DiracDelta(lambdaPrime, lambdaVal - 1) * DiracDelta(sigmaPrime, sigma + 1)
    )
    val = val * DiracDelta(lambdaPrime, lambdaVal)
    val = val + DiracDelta(sigmaPrime, sigma - 1)

    temp = (
        lambdaVal
        * sigma
        * DiracDelta(lambdaPrime, lambdaVal)
        * DiracDelta(sigma, sigmaPrime)
    )

    val = val1 + val + temp
    return val


#     ----- UTILS -----     #


def ConstantValues(hcutomega00: float, k: int, nu: float) -> tuple[float, float]:
    C: float = -2 * k * hcutomega00
    D: float = C * nu / 2
    return C, D


def DiracDelta(var1, var2) -> int:
    return 1 if var1 == var2 else 0


def SphericalHarmonics(l, lprime, lambdaVal, lambdaPrime) -> float:
    spherical = np.sqrt((5 * (2 * l + 1) / (4 * np.pi * 2 * lprime + 1)))
    spherical = spherical * cal(l, 2, [lambdaVal, 0], [lprime], [lambdaPrime])
    spherical = spherical * cal(l, 2, [0, 0], [lprime], [0])
    return spherical


def CalLambda(omega: float, sigma: list) -> list:
    lambdaVal: list = []
    for val in sigma:
        lambdaVal.append(omega - val)
    return lambdaVal


def CalLval(N: int) -> list:
    lVal = [n for n in range(0, N + 1, 2)]
    return lVal


def CalOmega(N: int) -> list:
    omegaMax = N + 0.5
    omega = []
    while omegaMax >= 0:
        omega.append(omegaMax)
        omegaMax = omegaMax - 1
    return omega


if __name__ == "__main__":
    sigma = [0.5, -0.5]
    main(sigma)
