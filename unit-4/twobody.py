from typing import List
import numpy as np
from CGcoeff import cal

# TODO Create a function to Calculate the T values


def main(
    j1: float, j2: float, j3: float, j4: float, n1: int, n2: int, n3: int, n4: int
) -> None:

    AT = 1
    jeve = True if j1 == j2 or j3 == j4 else False
    print(jeve)

    J1 = Jvalues(j1, j2)
    J2 = Jvalues(j3, j4)
    Jval = list(set(J1).intersection(set(J2)))
    M = Jvalues(1 / 2, 1 / 2)
    cgval = CGCal(j1, 1 / 2, j2, 1 / 2, J1, M)

    T = equivalent(jeve, Jval)

    for jval in Jval:
        print(jval)


def Calculate(
    j1: float,
    j2: float,
    j3: float,
    j4: float,
    n1: int,
    n2: int,
    n3: int,
    n4: int,
    J: float,
    AT: float,
) -> float:

    return j1


def SqrtElement(ja: float, jb: float, jc: float, jd: float) -> float:

    return np.sqrt(
        ((2 * ja + 1) * (2 * jb + 1) * (2 * jc + 1) * (2 * jd + 1))
        / ((1 + DiracDelta(ja, jb) * (1 + DiracDelta(jc, jd))))
    )


# ------UTILS ---------


def DiracDelta(Num1: float, Num2: float):
    return 1 if Num1 == Num2 else 0


def Jvalues(j1: float, j2: float) -> List:
    J = []
    Jmax = int(j1 + j2)
    Jmin = int(np.abs(j1 - j2))
    for i in range(Jmin, Jmax + 1, 1):
        J.append(i)

    return J


def CGCal(j1: float, m1: float, j2: float, m2: float, J: list, M: list) -> list:
    m1_m2 = [[m1, m2]]
    coeff: list = cal(j1, j2, m1_m2, M, J)
    return coeff


def equivalent(jeve: bool, J: list) -> list:
    T: list = []
    if jeve:
        for element, ind in enumerate(J):
            T.append(1 if J[ind] % 2 == 0 else 0)

    else:
        for element, ind in enumerate(J):
            T.append([0, 1])
    print(T)
    return T


if __name__ == "__main__":

    J1, J2 = 1 / 2, 3 / 2

    j1 = 3 / 2
    j2 = 3 / 2
    j3 = 3 / 2
    j4 = 3 / 2

    n1 = 0
    n2 = 0
    n3 = 0
    n4 = 0

    main(j1, j2, j3, j4, n1, n2, n3, n4)
