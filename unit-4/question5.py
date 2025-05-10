import math
import matplotlib.pyplot as plt
import numpy as np
from sympy import Matrix, factorial

def cgc(j1,j2,m1,m2,j):
    m=m1+m2
    if (j<=j1+j2 and j>=abs(j1-j2) and abs(m)<=abs(j) and abs(m1)<=abs(j1) and abs(m2)<=abs(j2)):
        k=0
        sum=0.0    
        cbcoeff=math.sqrt(((2*j+1)*factorial(int(j+j1-j2))*factorial(int(j-j1+j2))*factorial(int(j1+j2-j)))/factorial(int(j1+j2+j+1)))
        cbcoeff=cbcoeff*math.sqrt(factorial(int(j+m))*factorial(int(j-m))*factorial(int(j1-m1))*factorial(int(j1+m1))*factorial(int(j2-m2))*factorial(int(j2+m2)))
        for k in np.arange(max(-j+j2-m1,j1+m2-j,0),min(j1+j2-j,j1-m1,j2+m2)+1,1):
            sum=sum+pow((-1),k)*(pow(factorial(int(k))*factorial(int(j1+j2-j-k))*factorial(int(j1-m1-k))*factorial(int(j2+m2-k))*factorial(int(j-j2+m1+k))*factorial(int(j-j1-m2+k)),-1))
        
        return sum*cbcoeff
    else:
        return 0

def function1(A, delta_list, N_max):    
    Energy = []
    kap = 0.05
    mu = [0,0,0,0.35,0.625,0.63,0.448,0.434]
    for N in range (N_max+1):
        for i in range (N+1):
            om = i+0.5
            l_list = []
            lam_list = []
            sig_list = []
            H = []

            # Calculate basis for given value of N and om
            for l in range (N,-1,-2):
                for lam in range (-l, l+1):
                    sig = om - lam
                    if (abs(abs(sig)-0.5)<0.001):
                        l_list.append(l)
                        lam_list.append(lam)
                        sig_list.append(sig)

            # Calculate the hamiltonian matrix, <i|H|j>
            nbas = len(l_list)
            for delta in delta_list:
                H.append([])
                fdel = ((((1+2*delta/3)**(2))*(1-4*delta/3))**(-1/6))
                hw00 = 41*(A**(-1/3))
                hw0 = hw00*fdel
                C = -kap*2*hw00
                D = C*mu[N]/2.0
                for i in range(nbas):
                    H[-1].append([])
                    li = l_list[i]
                    lami = lam_list[i]
                    sigi = sig_list[i]
                    for j in range (nbas):
                        lj = l_list[j]
                        lamj = lam_list[j]
                        sigj = sig_list[j]

                        if (i==j):
                            h00 = (N+1.5)*hw0
                            hl2 = D*lj*(lj+1)
                            hls = C*lamj*sigj
                            hr2 = N+1.5
                        else:
                            h00 = 0
                            hl2 = 0
                        if (abs(li-lj)<0.001):
                            if (abs(lami-lamj-1)<0.001)and (abs(sigi-sigj+1)<0.001):
                                hls = C*0.5*math.sqrt((lj-lamj)*(lj+lamj+1))
                            elif (abs(lami-lamj+1)<0.001)and (abs(sigi-sigj-1)<0.001):
                                hls = C*0.5*math.sqrt((lj+lamj)*(lj-lamj+1))
                        else:
                            hls = 0
                        if (abs(lami-lamj)<0.001) and (abs(sigi-sigj)<0.001):
                            hY0 = cgc(lj, 2, lamj, 0, li)*cgc(lj, 2, 0, 0, li)*math.sqrt((2*lj+1)/(2*li+1))
                            if (abs(li-lj+2)<0.001):
                                hr2 = math.sqrt((N-lj+2)*(N+lj+1))
                            elif (abs(li-lj-2)<0.001):
                                hr2 = math.sqrt((N-lj)*(N+lj+3))
                        else:
                            hr2 = 0
                            hY0 = 0
                        hdelta = -delta*hw0*(2/3)*hr2*hY0
                        H[-1][-1].append(h00+hdelta+hls+hl2)
                        
            # Diagonalise H
            diagH = []
            diagHs = []
            for i in range (len(H)):
                P, diag = Matrix(H[i]).diagonalize()
                diagHs.append(sorted(np.diagonal(diag).tolist()))
            diagHs = np.array(diagHs)

            for i in range(len(l_list)):
                y = diagHs[:,i].tolist()
                Energy.append(y)
    return Energy

Z = 40
A = 84
delta_list = [0.2]
N_max = 7
Energy = function1(A, delta_list, N_max)
Energy = np.array(Energy)

N = A - Z
hw0 = 41*(A**(-1/3))

Energies = sorted(Energy[:,0])

# Defining Function
def f(x, T,Ener=Energies, Np=44):
    s = 0
    #flag = 1
    for i in range(len(Ener)):
        ek = Ener[i]
        dr = 1 + math.exp((ek-x)/T)
        s = s + 1/dr
    s = s - Np
    return s

# Defining derivative of function
def g(x, T,Ener=Energies, Np=40):
    s = 0
    #flag = 1
    for i in range(Np):
        ek = Ener[i]
        s = s+ 1/(4*T*((math.cosh((x-ek)/(2*T)))**2))
    return s

# Implementing Newton Raphson Method
def newtonRaphson(x0,T,N):
    step = 1
    flag = 1
    e=0.001
    eta=0.1
    condition = True
    while condition:
        x1 = x0 - eta*f(x0,T)/g(x0,T)
        x0 = x1
        step = step + 1

        if step > N:
            flag = 0
            break        
        
        condition = abs(f(x1,T)) > e
        if flag==1:
            return x1
        else:
            print('\n Not Convergent.')

# Fermi-Dirac including particle number

T_list = [0.5,1.0,1.5,2.0,2.5,3.0,3.5,4.0,4.5,5.0]
fig, (ax1, ax2) =  plt.subplots(1,2, figsize=(15, 6))
fig.suptitle('$^{84}Zr$ at $\delta = 0.2$ and T in MeV')

for T in T_list:
    n = []
    s=[]
    x0 = 44  # Starting guess
    Np = len(Energies)
    lam= newtonRaphson(x0, T,100)
    for i in range(Np):
        n.append(1/(1+math.exp((Energies[i]-lam)/T)))
        s.append(-( n[-1]*np.log(n[-1]) + (1-n[-1])*np.log(1-n[-1]) ))

    ax1.plot(Energies,n,label='{}'.format(T))
    ax2.plot(Energies,s,label='{}'.format(T))
    
ax1.legend(loc='upper right')
ax2.legend(loc='upper right')
ax1.set_ylabel('$n_{i}$')
ax2.set_ylabel('$s_{i}$')
ax1.set_xlabel('$e_{i}\; MeV$')
ax2.set_xlabel('$e_{i}\; MeV$')
#plt.axhline(0.5)
