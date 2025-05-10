import numpy as np
from scipy import integrate
from scipy.integrate import odeint, simpson
from scipy.optimize import bisect
import matplotlib.pyplot as plt

L, m, hbar, w, k, E0, En = 6, 1, 1, 1, 1, 0, 9
Erange = np.arange(E0, En, 0.1)
N = 1000
x, h = np.linspace(-L, L, N, retstep = True)
def V(x): return 0.5*k*x**2 # define potential   
psi = np.zeros(N)   

u = [0, 0.1]  # Initial conditions for wavefunction and its derivative

def f(u, x, E):
    y, z = u
    f1 = z
    f2 = ((2 * m) / (hbar**2)) * (V(x) - E) * y
    return (f1, f2)

def psiboundval(E):
    sol = odeint(f, u, x, args=(E,))
    return sol[:, 0][-1]  # Boundary condition at r=L

def shoot(Erange):
    Y = np.array([psiboundval(E) for E in Erange])
    eigval = np.array([bisect(psiboundval, Erange[i], Erange[i + 1])
                       for i in np.where(np.diff(np.signbit(Y)))[0]])
    return eigval

# Find the energy eigenvalues
En = shoot(Erange)
print("Energy Eigen Values:", En)

# Plot the normalized wavefunctions for each energy eigenvalue
#plt.figure(figsize=(10,15))
st = 0
# obtain the wave functions
for eigen in En:
    sol = odeint(f, u, x, args=(eigen,))
    psi = sol[:, 0]
    normpsi = psi/np.sqrt(simpson(psi*psi, x)) #Normalize the wave function
    #plt.subplot(int(len(En) / 2) + 1, 2, st+1)
    plt.ylabel(r"$\Psi_"+str(st)+"(x)$", fontsize=16)
    plt.xlabel(r'$x$', fontsize=16)
    plt.axhline(y = eigen, color="black")
    plt.axvline(color="black")
    plt.xlim(-L,L)
    plt.ylim(-1, 10)
    plt.plot(x, eigen + normpsi, label='$\Psi_'+str(st)+'(x)$') #plot wave function
    plt.plot(x, V(x))
    st = st + 1

plt.legend()
plt.tight_layout()
plt.show()
