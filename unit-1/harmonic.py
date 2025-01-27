from hermite import hermite
from factorial import factorial

import matplotlib.pyplot as plt
import numpy as np


def main() -> None:
    val = np.linspace(-8, 8, 100000)
    num = int(input("Enter the number for which you want to calculate\n"))
    plt.ylim(-3, 3)

    vharm = np.vectorize(harmonic)
    nrange = []

    for i in range(num):
        nrange.append(str(i))
        plt.plot(val, vharm(i, val))

    plt.legend(nrange)
    plt.show()


def harmonic(n: int, x: float) -> float:
    sig = (
        (1 / np.sqrt(2**n * factorial(n)))
        * np.power((22 / 7), (1 / 4))
        * np.exp(-(x**2 / 2))
        * hermite(n, x)
    )
    return sig


if __name__ == "__main__":
    main()
