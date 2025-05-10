"""

This is for sd shell

"""

import numpy as np
from sympy.physics.quantum.cg import CG
from sympy import S

def delta(a, b):
    """Kronecker delta function."""
    return 1.0 if a == b else 0.0

def get_n(l):
    """Returns the occupation number n based on the orbital angular momentum l."""
    return 1 if l == 2 else 2 if l == 0 else 1

def calculate_tbme(ja, jb, jc, jd, la, lb, lc, ld):
    """Calculates two-body matrix element (TBME) for SDI in the sd-shell with given conditions."""
    
    A_T = 1.0
    na, nb, nc, nd = get_n(la), get_n(lb), get_n(lc), get_n(ld)

    J1_range = range(int(abs(ja - jb)), int(ja + jb + 1))
    J2_range = range(int(abs(jc - jd)), int(jc + jd + 1))
    common_J_values = list(set(J1_range) & set(J2_range))

    results = {}

    for J in common_J_values:
        for T in [0, 1]:
            if (ja == jb or jc == jd) and (J + T) % 2 == 0:
                continue  # Skip if J+T is even for identical particles

            term1 = (-1)**(na + nb + nc + nd) * (A_T / (2 * (2 * J + 1))) * \
                    np.sqrt((2 * ja + 1) * (2 * jb + 1) * (2 * jc + 1) * (2 * jd + 1) / 
                            ((1 + delta(ja, jb)) * (1 + delta(jc, jd))))

            term2 = (-1)**(jb + jd + lb + ld) * \
                    CG(S(jb), S(-0.5), S(ja), S(0.5), J, 0).doit().evalf() * \
                    CG(S(jd), S(-0.5), S(jc), S(0.5), J, 0).doit().evalf() * \
                    (1 - (-1)**(la + lb + J + T))

            term3 = CG(S(jb), S(0.5), S(ja), S(0.5), J, 1).doit().evalf() * \
                    CG(S(jd), S(0.5), S(jc), S(0.5), J, 1).doit().evalf() * \
                    (1 + (-1)**(T))

            results[(J, T)] = round(term1 * (term2 - term3), 6)  # Round to 6 decimal places

    return results

# Define sd-shell quantum numbers
sd_shell = [(5/2, 2), (1/2, 0), (3/2, 2)]  # (j, l)

def compute_tbme(shell):
    """Computes and prints TBME values for the sd-shell."""
    print("\033[1mTBME for sd-shell:\033[0m")
    calculated_pairs = set()

    for ja, la in shell:
        for jb, lb in shell:
            for jc, lc in shell:
                for jd, ld in shell:
                    # Add conditions for interchanged pairs
                    pair1 = tuple(sorted([(ja, jb), (jc, jd)]))
                pair2 = tuple(sorted([(jc, jd), (ja, jb)]))
                pair3 = tuple(sorted([(jb, ja), (jc, jd)]))
                pair4 = tuple(sorted([(ja, jb), (jd, jc)]))
                pair5 = tuple(sorted([(jb, ja), (jd, jc)]))

                if (pair1 not in calculated_pairs and
                    pair2 not in calculated_pairs and
                    pair3 not in calculated_pairs and
                    pair4 not in calculated_pairs and
                    pair5 not in calculated_pairs):

                    results = calculate_tbme(ja, jb, jc, jd, la, lb, lc, ld)
                    print(f"<ja={ja},jb={jb}|V|jc={jc},jd={jd}>:")
                    for (J, T), result in results.items():
                        print(f"  J={J}, T={T}: {result:.6f}")
                    calculated_pairs.add(pair1)

# Compute TBME for sd-shell
compute_tbme(sd_shell)
