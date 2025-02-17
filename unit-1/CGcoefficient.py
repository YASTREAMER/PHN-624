import numpy as np
import pandas as pd

from factorial import factorial


def main() -> None:
    j1 = 0.5
    j2 = 0.5

    m1_m2 = [[0.5, 0.5], [-0.5, 0.5], [0.5, -0.5], [-0.5, -0.5]]

    M = [1, 0, 0, -1]
    J = [1, 1, 0, 1]

    cg = cal(j1, j2, m1_m2, M, J)

    print(f"The CGCoefficient is \n {pd.DataFrame(cg)}")


def cal(j1:float, j2:float, m1_m2:list, M:list, J:list) -> list:

    cg = []
    cg1 = []

    for i in range(len(M)):
        cg1 = []
        for j in range(len(m1_m2)):

            if m1_m2[j][0] + m1_m2[j][1] != M[i]:
                cg1.append(0.0)
                continue

            cg1.append(cgCoeff(j1, j2, J[i], m1_m2[j][0], m1_m2[j][1], M[i]))
        cg.append(cg1)

    return cg


def cgCoeff(j1, j2, J, m1, m2, M) -> tuple:

    F1 = np.sqrt(
        (
            factorial(j1 + j2 - J)
            * factorial(J + j1 - j2)
            * factorial(J + j2 - j1)
            * (2 * J + 1)
        )
        / factorial(j1 + J + j2 + 1)
    )

    F2 = np.sqrt(
        factorial(J + M)
        * factorial(J - M)
        * factorial(j1 + m1)
        * factorial(j1 - m1)
        * factorial(j2 + m2)
        * factorial(j2 - m2)
    )

    F3 = f3(j1, j2, J, m1, m2)

    return F1 * F2 * F3


def f3(j1, j2, J, m1, m2) -> float:

    a = j1 + j2 - J
    b = j2 + m2
    c = j1 - m1
    d = J - j2 + m1
    e = J - j1 - m2

    a, b, c, d, e = int(a), int(b), int(c), int(d), int(e)

    Smax = min(a, b, c)
    Smin = np.abs(min(d, e)) if min(d, e) < 0 else 0

    Smin, Smax = int(Smin), int(Smax)
    F3 = 0

    for s in range(Smin, Smax + 1):
        F3 = F3 + (
            (
                ((-1) ** s)
                * factorial(j1 - m1 - s)
                * factorial(j2 + m2 - s)
                * factorial(J - j2 + m1 + s)
            )
            * 1
            / ((factorial(J - j1 - m2 + s) * factorial(j1 + j2 - J - s) * factorial(s)))
        )
        # print(F3)
    return F3


if __name__ == "__main__":
    main()
