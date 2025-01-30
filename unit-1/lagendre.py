import numpy as np
import matplotlib.pyplot as plt
import os

from factorial import *


def main():
    dir()
    val = np.linspace(-1, 1, 100000)  # val is x
    num = int(
        input("enter the number for which you want to calculate\n")
    )  # num is also l

    m = []

    for i in range(-num, num):
        m.append(i)

    # plt.xlim(-2, 2)
    plt.ylim(-2, 2)

    nrange = []

    vlag = np.vectorize(associatedLagendre)

    for i in range(num):
        for j in m:
            plt.plot(val, vlag(i, j, val))
            nrange.append(str(i) + str(j))
            plt.legend(nrange)
        plt.savefig(f"Figure/associatedLagendre-{i}.png")
        plt.clf()
    # plt.show()


def lagendre(num: int, x: float) -> float:
    if num == 0:
        return 1
    elif num == 1:
        return x

    return (2 * num - 1) / num * (x * lagendre(num - 1, x)) - (
        (num - 1) / num * lagendre(num - 2, x)
    )


def associatedLagendre(l: int, m: int, x: float) -> float:

    if l <= 0:
        return 1

    if m < 0:
        return (
            -(1 ** np.exp(np.abs(m)))
            * (factorial(l - np.abs(m)) / factorial(l + np.abs(m)))
            * associatedLagendre(l, np.abs(m), x)
        )

    elif l == m:
        return (-(1**l)) * doublefact(2 * l - 1) * (1 - x**2) ** (l / 2)

    elif l == m + 1:
        return x * (2 * l - 1) * associatedLagendre(l - 1, l - 1, x)

    else:
        return (
            x * (2 * l - 1) * associatedLagendre(l - 1, m, x)
            - (l + m - 1) * associatedLagendre(l - 2, m, x)
        ) / (l - m)


def dir()->None:
    if not os.path.exists("Figure"):
        os.mkdir("Figure")

if __name__ == "__main__":
    main()
