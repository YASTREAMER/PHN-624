import matplotlib.pyplot as plt
import numpy as np


def main():

    potential = []
    rRange = []
    a = 0.5
    Vo = -25
    R = 5

    i = -10
    while not (i == 11):
        pot = wood(i, Vo, R, a)
        potential.append(pot)
        rRange.append(i)
        i += 1

    plt.plot(rRange, potential)
    plt.show()


def wood(r: int, Vo: int, R: int, a: float) -> float:
    r = np.abs(r)

    pot: float = Vo / (1 + (np.exp((r - R) / a)))

    return pot


if __name__ == "__main__":
    main()
