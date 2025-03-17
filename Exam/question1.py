import numpy as np
import os



def main():
    k = int(input("Enter the value of k \n"))
    l = int(input("Enter the value of l \n"))
    div = 10000
    result = integration(harmoinc, 100, -100, div, k, l)
    result = np.round(result)
    result = np.abs(result)

    if os.path.isfile("Question1.txt"):
        os.remove("Question1.txt")

    file = open("Question1.txt", "a")
    file.write(f"The value of the integral if {result} \n")
    print(f"The value of the integral for value of {k} and {l} is {result} \n")


def harmoinc(k, l, x) -> float:
    result = wavefunction(k, x) * wavefunction(l, x)
    return result


def wavefunction(n: int, x: float) -> float:

    wavefunc = 1 / np.sqrt(2**n * factorial(n))
    wavefunc = wavefunc / ((np.pi) ** (1 / 4))
    wavefunc = wavefunc * np.exp((-(x**2)) / 2) * hermite(n, x)
    return wavefunc


def factorial(n: int) -> int:
    if n == 0 or n == 1:
        return 1

    if n < 0:
        return -factorial(np.abs(n))

    return n * factorial(n - 1)


def hermite(num: int, x: float) -> float:

    if num == 0:
        return 1
    elif num == 1:
        return 2 * x
    return 2 * x * hermite(num - 1, x) - 2 * (num - 1) * hermite(num - 2, x)


def integration(fx, upper, lower, div=1000, *args) -> float:

    result = 0
    num = 1

    for i in np.linspace(lower, upper, div):

        if i == lower or i == upper:
            result = result + fx(*args, i)
            num += 1
            continue

        elif num % 2 == 0:
            result = result + 2 * fx(*args, i)

        else:
            result = result + 4 * fx(*args, i)

        num += 1

    x = (upper - lower) / div
    result = (x / 3) * result

    return result


if __name__ == "__main__":
    main()
