from typing import List
import numpy as np

# TODO Create a function to calculate the Jvalues


def main():
    j1 = 5 / 2
    j2 = 5 / 2
    j3 = 5 / 2
    j4 = 5 / 2

    J1 = Jvalues(j1, j2)
    print(J1)


def matrix(ja: float, jb: float, jc: float, jd: float) -> float:

    return 0.1


def JvalueDef(J1: list, J2: list) -> list:
    Jval = list(set(J1).intersection(set(J2)))
    return Jval


def Jvalues(j1: float, j2: float) -> List:
    J = []
    Jmax = int(j1 + j2)
    Jmin = int(np.abs(j1 - j2))
    for i in range(Jmin, Jmax + 1, 1):
        J.append(i)

    return J


if __name__ == "__main__":
    main()
