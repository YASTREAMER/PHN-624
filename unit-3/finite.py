import numpy as np
import matplotlib.pyplot as plt

from util import runge


def main():
    lower = -12
    upper = 12
    Z = [0, 0.1]

    temp = runge(Z, potential, lower, upper)
    plt.plot(temp[0], temp[1])
    plt.show()


def potential(x: int):

    if np.abs(x) <= 10:
        return 25
    else:
        return 0


if __name__ == "__main__":
    main()
