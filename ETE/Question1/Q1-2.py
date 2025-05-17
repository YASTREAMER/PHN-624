import numpy as np
import matplotlib.pyplot as plt

h_cut = 197.327  # in MeV-fm
m = 938.272  # in MeV/c^2
V_0 = -38  # in MeV
L = 10  # in fm
h = (L - (-L)) / 100
x = np.arange(-L, L, h)
Z1 = np.empty(np.size(x))
Z2 = np.empty(np.size(x))
a = 0.5
R = 6

Z1[0] = 0
Z2[0] = 1.0


def WoodSaxon(d):
    return V_0 / (1 + np.exp(abs(d) - R) / a)


def func1(x, z1, z2):
    return z2


def func2(x, z1, z2, V, E):
    return (2 * m / h_cut**2) * (V - E) * z1


def rk4(x, Z1, Z2, function, e):
    for j in range(np.size(x) - 1):
        k_11 = h * func1(x[j], Z1[j], Z2[j])
        k_21 = h * func2(x[j], Z1[j], Z2[j], function(x[j]), e)

        k_12 = h * func1(x[j] + h / 2, Z1[j] + k_11 / 2, Z2[j] + k_21 / 2)
        k_22 = h * func2(
            x[j] + h / 2, Z1[j] + k_11 / 2, Z2[j] + k_21 / 2, function(x[j] + h / 2), e
        )

        k_13 = h * func1(x[j] + h / 2, Z1[j] + k_12 / 2, Z2[j] + k_22 / 2)
        k_23 = h * func2(
            x[j] + h / 2, Z1[j] + k_12 / 2, Z2[j] + k_22 / 2, function(x[j] + h / 2), e
        )

        k_14 = h * func1(x[j] + h, Z1[j] + k_13, Z2[j] + k_23)
        k_24 = h * func2(x[j] + h, Z1[j] + k_13, Z2[j] + k_23, function(x[j] + h), e)

        Z1[j + 1] = Z1[j] + (k_11 + 2 * (k_12 + k_13) + k_14) / 6
        Z2[j + 1] = Z2[j] + (k_21 + 2 * (k_22 + k_23) + k_24) / 6
    return Z1


def f(e):
    z = rk4(x, Z1, Z2, WoodSaxon, e)
    return z[-1]


def bisection(E_low, E_high, tol=1e-6, max_iter=100):
    if f(E_low) * f(E_high) > 0:
        return None
    iter_count = 0
    while abs(E_high - E_low) > tol and iter_count < max_iter:
        E_mid = (E_low + E_high) / 2
        f_mid = f(E_mid)

        if abs(f_mid) < tol:
            return E_mid

        if f(E_low) * f_mid < 0:
            E_high = E_mid
        else:
            E_low = E_mid

        iter_count += 1

    return (E_low + E_high) / 2


def find_eigenvalues(E_min, E_max, num_states, tol=1e-6):
    eigenvalues = []
    E_low = E_min
    step = (E_max - E_min) / (10000)

    while len(eigenvalues) < num_states and E_low < E_max:
        E_high = E_low + step

        if f(E_low) * f(E_high) < 0:
            eigen_energy = bisection(E_low, E_high, tol)
            if eigen_energy is not None:
                eigenvalues.append(eigen_energy)

        E_low = E_high

    return eigenvalues


E_min = -50
E_max = 50
num_states = 4
eigenvalues = find_eigenvalues(E_min, E_max, num_states)

for en in eigenvalues:
    Z = rk4(x, Z1, Z2, WoodSaxon, en)
    plt.plot(x, Z)
plt.title("Eigen states in WS potential")
plt.grid()
plt.savefig("Wood-Saxon.png")
plt.show()
