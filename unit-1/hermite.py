import matplotlib.pyplot as plt
import numpy as np


def hermite(num: int, x: float) -> float:
    if num == 0:
        return 1
    elif num == 1:
        return 2 * x
    return 2 * x * hermite(num - 1, x) - 2 * (num - 1) * hermite(num - 2, x)


def main():
    val = np.linspace(-3, 3, 100000)
    num = int(input("enter the number for which you want to calculate\n"))
    plt.ylim(-50, 50)

    nrange=[]

    vherm = np.vectorize(hermite)

    for i in range(num):
        plt.plot(val, vherm(i, val))
        nrange.append(str(i))

    plt.legend(nrange)
    plt.show()


if __name__ == "__main__":
    main()
