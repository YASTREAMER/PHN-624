import numpy as np


def main():
    A = 36
    hw = 41 * A ** (-1 / 3)
    oscLen = osclen(hw)
    hW = 45 * A ** (-1 / 3) - 25 * A ** (-2 / 3)
    OscLen = osclen(hW)
    print(f"The first oscillator lenght is {oscLen}")
    print(f"The second oscillator lenght is {OscLen}")
    dif = OscLen - oscLen
    print(f"The difference between them is {dif}")


def osclen(hw: float) -> float:
    oscLen = 197.33 / np.sqrt(940 * hw)
    return oscLen


if __name__ == "__main__":
    main()
