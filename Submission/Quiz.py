"""

This is the code for the quiz of quadrapole moment and magentic moment.

Since there will never be contribution from gnL thus we can ignore that term

"""


def main() -> None:

    # Constants
    j = 3 / 2
    gl = 1
    gs = 5.585

    num = int(input("Enter a atomic number for Lithium \n"))

    radius = rad(num)
    quad = quadrapole(radius, j)
    print(f"The quadrapole moment is {quad}")

    magmo = magneticMoment(j, gl, gs)
    print(f"The magnetic moment is {magmo}")


def rad(num: float) -> float:
    radius = (3 / 5) * ((1.2 * num ** (1 / 3)) ** 2)
    return radius


def quadrapole(radius: float, j: float) -> float:

    quad = -((2 * j - 1) / (2 * (j + 1))) * radius
    return quad


def magneticMoment(j: float, gl: float, gs: float) -> float:

    magMoment = (j - 0.5) * gl + 0.5 * gs

    return magMoment


if __name__ == "__main__":
    main()
