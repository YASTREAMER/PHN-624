import numpy as np
from spherical import *


def main():

    theta = np.linspace(-np.pi, np.pi, 1000)  # Adjusting the resolution of theta
    phi = np.linspace(0, 2 * np.pi, 1000)

    theta, phi = np.meshgrid(theta, phi)  # Create mesh grid for theta and phi
    m = 2
    l = 5
    temp = Cal(l, m, theta, phi)
    print(temp)


def hydrogenCal(n: int, l: int, m: int, theta: float, phi: float):

    spherical = Cal(l, m, theta, phi)
    return spherical


if __name__ == "__main__":
    main()
