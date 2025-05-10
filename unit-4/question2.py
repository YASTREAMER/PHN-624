import numpy as np
import matplotlib.pyplot as plt

def energy_eigen_value(N,delta,n_z):
    return (N + 3/2 - (1/3) * delta * (3 * n_z - N))

N_max = int(input("Enter the maximum value of N : "))

delta = np.linspace(-1,1,100)
N = np.arange(0,N_max+1,1)
#n_z = np.arange(0,N+1,1)

for n in N:
    for i in range(n):
        E = [energy_eigen_value(n,d,i) for d in delta]
        plt.plot(delta,E,label=f'N={n:.2f}')
plt.title("Anisotropic harmonic oscillator")
plt.xlabel("δ")
plt.ylabel("E")
# plt.ylim(1,7)
# plt.xlim(-1,1)
# plt.legend()
plt.grid()
plt.show()

A = 80

def new_energy_eigen_value(A, N, delta, n_z):
    factor1 = 41 * A**(-1/3)
    factor2 = ((1 + (2/3)*delta) **2 * (1 - (4/3)*delta)) **(-1/6)
    factor3 = (N + 3/2 - (1/3) * delta * (3*n_z - N))
    return factor1 * factor2 * factor3

for n in N:
    for i in range(n):
        E_new = [new_energy_eigen_value(A,n,d,i) for d in delta]
        plt.plot(delta,E_new,label=f'N={n:.2f}')
plt.xlabel("δ")
plt.ylabel("E")
# plt.ylim(10,70)
# plt.xlim(-0.8,0.8)
# # plt.legend()
plt.grid()
plt.show()
