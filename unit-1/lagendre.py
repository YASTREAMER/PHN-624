import numpy as np
import matplotlib.pyplot as plt

from factorial import *


def main():
    plt.ylim(-1, 1)
    val = np.linspace(-1, 1, 100000)  # val is x
    choice = int(input("Enter 1 for Legendre and 2 for Associated Legendre"))
    num = int(
        input("enter the number for which you want to calculate\n")
    )  # num is also l

    if choice == 2:
        m = num
        plotAssocaite(num, val, m)
    elif choice == 1:
        plotLeg(num, val)
    else:
        print("Error: Invalid input")


def plotLeg(num, val) -> None:

    nrange = []
    vlag = np.vectorize(lagendre)

    for i in range(num):
        plt.plot(val, vlag(i, val))
        nrange.append(str(i))
    plt.legend(nrange)
    plt.title("Legendre")
    plt.show()


def plotAssocaite(num: int, val, m: int) -> None:

    nrange = []
    vlag = np.vectorize(associatedLagendre)

    for i in range(m + 1):
        plt.plot(val, vlag(num, i, val))
        nrange.append(str(num) + str(i))
    plt.legend(nrange)
    plt.title("Associated Legendre")
    plt.show()


def lagendre(num: int, x: float) -> float:

    if num == 0:
        return 1
    elif num == 1:
        return x

    return (2 * num - 1) / num * (x * lagendre(num - 1, x)) - (
        (num - 1) / num * lagendre(num - 2, x)
    )


def associatedLagendre(l: int, m: int, x: float) -> float:

    if l == 0:
        return 1

    if m < 0:
        return (
            -(1 ** np.abs(m))
            * (factorial(l - np.abs(m)) / factorial(l + np.abs(m)))
            * associatedLagendre(l, np.abs(m), x)
        )

    elif l == m:
        return (-1) ** l * doublefact(2 * l - 1) * np.power(1 - x**2, l / 2)

    elif l == m + 1:
        return x * (2 * l - 1) * associatedLagendre(l - 1, l - 1, x)

    else:
        P = x * (2 * l - 1) * associatedLagendre(l - 1, m, x) - (
            l + m - 1
        ) * associatedLagendre(l - 2, m, x)
        return P / (l - m)


if __name__ == "__main__":
    main()
