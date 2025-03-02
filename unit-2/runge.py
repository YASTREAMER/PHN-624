import numpy as np
import matplotlib.pyplot as plt


def main():
    """
    The first step in runge kutta method is defining the initial condition. Since this code will also be used in unit-3, a template of
    code is very helpful

    """

    runge()


def runge(A=197, atomic=79):

    E = 5  # Mev
    # For an alpha particle
    m = 2 * 939.5654133 + 2 * 938.272083
    # m = m *1.6e-16

    RT = 1.2 * A ** (1 / 3)
    # Q = atomic * (1e-6)
    # q = 2e-6
    Q = 79
    q = 2

    """
     Defining the initial condition here
    """

    initx = -100  # fm
    initxdot = np.sqrt(2 * E / m)
    inity = -200  # fm
    initydot = 0  # fm
    plt.xlim(-400, 400)
    plt.ylim(-400, 400)

    for i in range(inity, inity + 401, 5):
        print(f"The iteration number is {i}")
        Z = [initx, initxdot, i, initydot]
        calculate(Z, q, Q, m, RT)

    plt.show()


def calculate(Z, q, Q, m, RT)-> None:

    t = 0
    step = 10000
    x = []
    y = []
    k1 = []
    k2 = []
    k3 = []
    k4 = []

    for i in range(0, step):
        k1 = []
        k2 = []
        k3 = []
        k4 = []

        rR = np.sqrt(Z[0] ** 2 + Z[2] ** 2)
        r = max(RT, rR)
        temp = Z

        for j in range(0, 4, 1):
            k1.append(f1(t, temp))
            k2.append(f2(t, temp, q, Q, m, r))
            k3.append(f3(t, temp))
            k4.append(f4(t, temp, q, Q, m, r))
            temp = []
            temp = [
                Z[0] + k1[j] / 2,
                Z[1] + k2[j] / 2,
                Z[2] + k3[j] / 2,
                Z[3] + k4[j] / 2,
            ]

        Z[0] = Z[0] + ((2 * (k1[1] + k1[2]) + k1[0] + k1[3]) / 6)

        Z[1] = Z[1] + ((2 * (k2[1] + k2[2]) + k2[0] + k2[3]) / 6)

        Z[2] = Z[2] + ((2 * (k3[1] + k3[2]) + k3[0] + k3[3]) / 6)

        Z[3] = Z[3] + ((2 * (k4[1] + k4[2]) + k4[0] + k4[3]) / 6)
        x.append(Z[0])
        y.append(Z[2])
        t = t + i
    plt.plot(x, y)


def f1(t: int, z: list) -> float:

    return z[1]


def f2(t: int, z: list, q, Q, m, r) -> float:

    return q * Q * 1.44 * z[0] / (m * r**3)


def f3(t: int, z: list) -> float:

    return z[3]


def f4(t: int, z: list, q, Q, m, r) -> float:

    return q * Q * 1.44 * z[2] / (m * r**3)


if __name__ == "__main__":
    main()
