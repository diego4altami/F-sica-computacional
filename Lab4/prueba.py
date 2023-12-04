import numpy as np
import numpy.random as nr
import matplotlib.pyplot as plt

NTl = 1000 #Numero de atomos de Tl
NPb = 0 #Numero de atomos de Pb
h = 1  # intervalo de tiempo
tau = 3.053*60  #tiempo de vida media
p = 1 - 2**(-h/tau) #probabilidad de decaimiento

NTalios = [NTl]
NPlomos = [NPb]

tmax = 1000 

for _ in range(1,tmax,h):

    decays = 0 #Cuenta el número de átomos que decayeron en el tiempo t
    for i in range(NTl):
        r = nr.random()
        if r < p: #Moneda sesgada
            decays += 1
    NTl -= decays
    NPb += decays
    NTalios.append(NTl)
    NPlomos.append(NPb)

tpoints = np.arange(0,tmax,h)
logNtalio=np.log(NTalios)

def funcion1(NTl,t,tau):
    Ntalfunc=NTl*2**(-t/tau)
    return Ntalfunc

def minimoscuadrados(x, y):
    X = np.vstack([x, np.ones(len(x))]).T
    m,b= np.linalg.lstsq(X, y, rcond=None)[0]
    return m,b

NTalios_func = funcion1(NTl, tpoints, tau) # Use a different variable name

pendiente,ordenada = minimoscuadrados(tpoints, NTalios)

NTalios_arr = np.array(NTalios)
fit = pendiente * NTalios_arr + ordenada

plt.plot(tpoints, NTalios, label='data')
plt.plot(tpoints, fit, label='fit')
plt.legend()
plt.show()
