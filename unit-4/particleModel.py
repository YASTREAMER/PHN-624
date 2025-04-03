import matplotlib.pyplot as plt


def main():
    Nmax = int(input("Enter the value for N\n"))
    for N in range(0, Nmax + 1):
        for l in range(N, -1, -2):
            en, en1 = model(N, l)
            plt.plot(en)
            plt.plot(en1)
            print(f"Energy val is {en}")
            print(f"Energy val is {en1}")
    plt.show()


def model(N: int, l: int) -> tuple:
    energy1 = []
    energy2 = []
    ldots = ldotsfunc(l)
    lsqaure = lsquarefunc(l)
    en = N + 1.5

    for i in range(0, 17, 1):
        if i < 9:
            energy1.append(en)
            energy2.append(en)
        else:
            energy1.append(en + lsqaure + ldots[0])
            energy2.append(en + lsqaure + ldots[1])

    return energy1, energy2


def ldotsfunc(l: int) -> list:
    ldots = []
    if l == 0:
        ldots.append(0)
        ldots.append(0)
        return ldots

    ldots.append(-0.1 * l / 2)
    ldots.append(0.1 * (l + 1) / 2)
    return ldots


def lsquarefunc(l: int) -> float:
    return l * (l + 1) * (-0.0225)


if __name__ == "__main__":
    main()
