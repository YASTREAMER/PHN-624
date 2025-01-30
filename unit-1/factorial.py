import numpy as np


def factorial(n: int) -> int:
    if n == 0 or n == 1:
        return 1

    if n < 0:
        return -factorial(np.abs(n))

    return n * factorial(n - 1)


def doublefact(n: int) -> int:
    if n == 0 or n == 1:
        return 1

    return n * doublefact(n - 2)


def main() -> None:

    n = int(input("Enter the number \n"))

    print(f"The factorial of {n} is {factorial(n)}")
    print(f"The double factorial of {n} is {doublefact(n)}")


if __name__ == "__main__":
    main()
