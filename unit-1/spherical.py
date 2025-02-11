import numpy as np

from lagendre import *
from factorial import *


def main():
    l = 0
    m = 0
    theta = np.linspace(-np.pi, np.pi, 1000)  # Adjusting the resolution of theta
    phi = np.linspace(0, 2 * np.pi, 1000)

    theta, phi = np.meshgrid(theta, phi)  # Create mesh grid for theta and phi

    x, y, z = Cal(l, m, theta, phi)

    fig = plt.figure(figsize=(8, 8))
    ax = fig.add_subplot(111, projection="3d")
    ax.set_xlabel("X")
    ax.set_ylabel("Y")
    ax.set_zlabel("Z")
    ax.plot_surface(x, y, z, cmap="viridis")
    plt.title("Real part")

    plt.show()


def Cal(l, m, theta, phi):
    r = Spherical(l, m, theta, phi)
    return (
        r * np.sin(theta) * np.cos(phi),
        r * np.sin(theta) * np.sin(phi),
        r * np.cos(theta),
    )


def Spherical(l: int, m: int, theta: float, phi: float) -> float:
    fp =np.sqrt((2 * l + 1) / (4 * np.pi) * (factorial(l - m) / factorial(l + m)))
    sphere = fp * associatedLagendre(l, m, np.cos(theta)) * np.cos(phi * m)
    return sphere**2


if __name__ == "__main__":
    main()
