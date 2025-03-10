import numpy as np
from CGcoefficient import cal


def main():
    j1, j2, m, J, M = (17 / 2), 2, [[2.5, 0]], [(13 / 2)], [2.5]
    trans = transition(j1, j2, m, M, J)
    print(f"The transition strength is for {j1} is  { trans }")

    j1, j2, m, J, M = (21 / 2), 2, [[2.5, 0]], [(17 / 2)], [2.5]
    trans = transition(j1, j2, m, M, J)
    print(f"The transition strength is for {j2} is  { trans }")


def transition(j1, j2, m, M, J):

    qt = 224

    cgcoeff = cal(j1, j2, m, M, J)

    w3 = (-1) ** (j1 - j2 - M[0]) / (2 * J[0] + 1) ** (1 / 2)

    temp = w3 * cgcoeff[0][0]

    trans = np.power(temp, 2) * qt * 5 / (16 * np.pi)
    print(trans)
    trans = trans * 1.6e-19 *1e-15

    return trans


if __name__ == "__main__":
    main()
