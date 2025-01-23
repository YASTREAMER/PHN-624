import numpy as np


def main():

    i = 0
    lr = 0.0001
    lowerlimit = 0
    upperlimit = 1
    integration = 0
    while i <= 2:
        hel = fun(i)
        if i == 2 or i == 0:
            integration = integration + hel
        else:
            integration = integration + (hel * 2)
        i += lr

    integration = (upperlimit - lowerlimit) * lr * integration * 0.5

    print(integration)


def fun(x: float) -> float:

    return x**2 * np.exp(-(x**2))


if __name__ == "__main__":
    main()
