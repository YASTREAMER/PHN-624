import numpy as np
import matplotlib.pyplot as plt

from lagendre import *
from factorial import *


def main():
    l = 3
    m = 2
    theta = np.linspace(0, np.pi, 100)  # Adjusting the resolution of theta
    phi = np.linspace(0, 2 * np.pi, 100)  # Adjusting the resolution of phi

    theta, phi = np.meshgrid(theta, phi)  # Create mesh grid for theta and phi

    x, y, z = realCal(l, m, theta, phi)

    fig = plt.figure()
    ax = fig.add_subplot(111, projection="3d")
    ax.set_xlabel("X")
    ax.set_ylabel("Y")
    ax.set_zlabel("Z")
    ax.plot_surface(x, y, z, cmap="viridis")
    plt.title("Real part")

    x, y, z = imagCal(l, m, theta, phi)

    fig = plt.figure()
    ax = fig.add_subplot(111, projection="3d")
    ax.set_xlabel("X")
    ax.set_ylabel("Y")
    ax.set_zlabel("Z")
    ax.plot_surface(x, y, z, cmap="viridis")
    plt.title("Imaginary part")
    

    plt.show()


def realX(theta, phi):
    return np.sin(theta) * np.cos(phi)


def realY(theta, phi):
    return np.sin(theta) * np.sin(phi)


def realZ(m, phi):
    return np.cos(np.abs(m) * phi)


def realCal(l, m, theta, phi):
    return (
        realX(theta, phi) * realSpherical(l, m, theta, phi),
        realY(theta, phi) * realSpherical(l, m, theta, phi),
        realZ(m, phi) * realSpherical(l, m, theta, phi),
    )

def imagCal(l, m, theta, phi):
    return (
        realX(theta, phi) * imagSpherical(l, m, theta, phi),
        realY(theta, phi) * imagSpherical(l, m, theta, phi),
        realZ(m, phi) * imagSpherical(l, m, theta, phi),
    )


def realSpherical(l: int, m: int, theta: float, phi: float) -> float:
    sphere = (
        np.sqrt(
            ((2 * l + 1) / 4 * (22 / 7))
            * (factorial(l - abs(m)) / factorial(l + abs(m)))
        )
        * associatedLagendre(l, m, np.cos(theta))
        * np.cos(m * phi)
    )
    return np.power(sphere, 2)

def imagSpherical(l: int, m: int, theta: float, phi: float) -> float:
    sphere = (
        np.sqrt(
            ((2 * l + 1) / 4 * (22 / 7))
            * (factorial(l - abs(m)) / factorial(l + abs(m)))
        )
        * associatedLagendre(l, m, np.cos(theta))
        * np.sin(m * phi)
    )
    return np.power(sphere, 2)


if __name__ == "__main__":
    main()
