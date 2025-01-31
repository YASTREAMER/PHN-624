import numpy as np

from lagendre import *
from factorial import *


def main():

    pass


def spherical(l: int, m: int, theta:int) -> float:
    return (
        np.sqrt(((2 * l - 1) / 2 * np.pi) * factorial(l - m) / factorial(l + m))
        * associatedLagendre(l,m,np.cos(theta)) * np.exp(m)
    )


if __name__ == "__main__":
    main()
