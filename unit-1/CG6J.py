from factorial import factorial

import numpy as np

def main():
    j1 = 2
    j2 = 0
    j3 = 0
    J1 = 1
    J2 = 0
    J3 = 0

    symbol = weigner(j1, j2, j3, J1, J2, J3)
    print(symbol)


def weigner(j1, j2, j3, J1, J2, J3) -> float:

    symbol = np.sqrt(
        triangle(j1, j2, j3)
        * triangle(j1, J2, J3)
        * triangle(J1, j2, J3)
        * triangle(J1, J2, j3)
        * funk(j1, j2, j3, J1, J2, J3)
    )

    return symbol


def funk(j1, j2, j3, J1, J2, J3) -> float:
    val = 0

    tmax = min(j1 + j2 + j3, j1 + J2 + J3, J1 + j2 + J3, J1 + J2 + j3)
    tmin = max(j1 + j2 + J1 + J2, j2 + j3 + J2 + J3, j1 + j3 + J1 + J3)

    for t in range(tmin, tmax + 1):

        val = (
            factorial(t - j1 - j2 - j3)
            * factorial(t - j1 - J2 - J3)
            * factorial(t - J1 - j2 - J3)
            * factorial(t - J1 - J2 - j3)
            * factorial(j1 + j2 + J1 + J2 - t)
            * factorial(j2 + j3 + J2 + J3 - t)
            * factorial(j1 + j3 + J1 + J3 - t)
        )

        val = ((-1) ** t) * factorial(t + 1) * val

    return val


def triangle(a, b, c) -> float:
    delta = factorial(a + b - c) * factorial(a - b + c) * factorial(-a + b + c)
    delta = delta / factorial(a + b + c + 1)
    return delta


if __name__ == "__main__":
    main()
