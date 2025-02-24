import numpy as np


def main():
    lowerlimit = 0
    upperlimit = 1
    integration(fun, upperlimit, lowerlimit)


def fun(x: float) -> float:

    return x**2 * np.exp(-(x**2))


def integration(fx, upper, lower, div=1000, *args) -> float:

    result = 0
    num = 1

    for i in np.linspace(lower, upper , div):

        if i == lower or i ==upper:
            print(i)
            result = result + fx(i, *args)
            num += 1
            continue

        elif num % 2 == 0:
            result = result + 2 * fx(i, *args)

        else:
            result = result + 4 * fx(i, *args)

        num += 1

    x = (upper - lower) / div

    result = (x / 3) * result

    return result


if __name__ == "__main__":
    main()
