def runge(Z: list, potentiaL, lower, upper) -> list:
    hcut = 197.3269631  # in MeV
    mass = 938.272  # in MeV
    k1 = []
    k2 = []

    ans = Z
    result = []

    for E in range(1, 20):
        Z = ans

        for x in range(lower, upper):
            V = potentiaL(x)
            k1 = []
            k2 = []

            temp = Z

            for j in range(0, 4, 1):
                k1.append(f1(temp))
                k2.append(f2(temp, mass, hcut, V, E))
                temp = []
                if j == 0 or j == 3:
                    temp = [
                        Z[0] + k1[j],
                        Z[1] + k2[j],
                    ]
                else:
                    temp = [
                        Z[0] + k1[j] / 2,
                        Z[1] + k2[j] / 2,
                    ]

            Z[0] = Z[0] + ((2 * (k1[1] + k1[2]) + k1[0] + k1[3]) / 6)

            Z[1] = Z[1] + ((2 * (k2[1] + k2[2]) + k2[0] + k2[3]) / 6)

        result.append(Z)

    return result


def f1(Z) -> float:
    return Z[1]


def f2(Z, m, hcut, V, E) -> float:
    val = 2 * m / hcut**2
    val = val * (V - E) * Z[1]
    return val
