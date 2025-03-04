import numpy as np
from CGcoefficient import cal


def main():
    j1,j2,m,J,M = (13/2), 2, [[2.5,2.5]],[(17/2)] , 2.5

    trans = transition(j1,j2,m,J,M)
    print(trans)


def transition(j1, j2, m, M, J):

    qt = 224

    cgcoeff = cal(j1, j2, m, M, J)

    trans = np.power(cgcoeff,2) * qt * 5 / (16 * np.pi)

    return trans

if __name__ == "__main__":
    main()
