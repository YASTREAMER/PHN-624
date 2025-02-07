import numpy as np
import matplotlib.pyplot as plt

# import os


def main():
    # dir()
    plt.ylim(-10, 10)
    val = np.linspace(-10, 10, 100000)  # val is x

    num = int(
        input("Enter the number for which you want to calculate\n")
    )  # num is also l

    nrange = []
    vlag = np.vectorize(laguerre)

    for i in range(num):
        plt.plot(val, vlag(i, val))
        nrange.append(str(i))
    plt.legend(nrange)
    plt.title("Laguerre")
    plt.show()

    valag = np.vectorize(associatedLaguerre)

    for i in range(num):

        if i > num:
            break

        plt.plot(val, valag(num, i, val))
        nrange.append(f"{str(num) + str(i)}")
    plt.legend(nrange)
    plt.title("Associated Laguerre")
    plt.show()


def laguerre(n: int, x: float) -> float:

    if n == 0:
        return 1
    elif n == 1:
        return -x + 1

    return ((2 * n - 1 - x) * laguerre(n - 1, x) - (n - 2) * laguerre(n - 2, x)) / n


def associatedLaguerre(n: int, k: int, x: float):

    if n == 0:
        return 1

    elif n == 1:
        return -x + k + 1

    return ((2 * n + k + 1 - x) / (n + 1)) * associatedLaguerre(n - 1, k, x) - (
        (n + k) / (n + 1)
    ) * associatedLaguerre(n - 2, k, x)


if __name__ == "__main__":
    main()
