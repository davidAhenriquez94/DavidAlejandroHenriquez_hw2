import numpy as np
import matplotlib.pyplot as plt

c = 299792458 
def planck(v,T):
    h = 6.62606957*(10**(-34))
    k = 1.3806488*(10**(-23))
    return (2*h*(v**3))/((c**2)*(np.exp((h*v)/(k*T)) -1))

def derivada_planck(v,T):
    dv = 10**(10)
    return (1.0/dv)*(planck(v+dv,T) - planck(v,T))

def segunda_derivada_planck(v,T):
    dv = 10**(13)
    return (1.0/(dv**2))*(planck(v+dv,T) + planck(v-dv,T) - 2*planck(v,T))

def maximo(T,x_0):
    x = x_0
    error = 0.001*(10**14)
    A = []
    r = 10**14
    while abs(r) > error:
        A.append(x)
        x_n = x - (derivada_planck(x,T)/segunda_derivada_planck(x,T))
        r = x_n - x
        x = x_n
    return A[-1]


s = []

for i in range(1,17):
    s.append(c/(maximo(5*10**(i/2.0),0.3*(10**(12+i/2.0))))*(5*10**(i/2.0)))
constante = str(np.mean(s))
fileout = open('b.dat','w')
fileout.write(constante)
fileout.close()
