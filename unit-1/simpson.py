import re
import numpy as np


def integration(fx, upper, lower, div=1000, *args) -> float:

    result = 0
    num = 1

    for i in np.linspace(lower, upper, div):

        if i == lower or i == upper:
            result = result + fx(*args, i)
            num += 1
            continue

        elif num % 2 == 0:
            result = result + 2 * fx(*args,i)

        else:
            result = result + 4 * fx(*args,i)

        num += 1

    x = (upper - lower) / div
    result = (x / 3) * result


    return result
