import numpy as np
import matplotlib.pyplot as plt

h = 6.626*(10**(-34))
c = 299792458.0
k = 1.380*(10**(-23)) 

def emision_de_cuerpo_negro(v,T):
   return ((2*h*(v**3))/(c**2))*(1.0/(np.exp((h*v)/(k*T)) - 1))
    
def integral(T,s_0):
    suma_impar = 0.0
    suma_par = 0.0
    s = s_0
    for i in np.arange(2,199):
       if i%2 == 1:
          suma_impar = suma_impar + emision_de_cuerpo_negro(i*s,T)
       else:
          suma_par = suma_par + emision_de_cuerpo_negro(s*i,T)
    integr = (1/3.0)*s*(emision_de_cuerpo_negro(s,T) + emision_de_cuerpo_negro(200*s,T) + 2*suma_impar + 4*suma_par)
    return integr

I = [] 
for i in range(2,17):
   I.append((integral(5*(10**(i/2.0)),0.01*(10**(13+i/2.0))))*(np.pi/((5*(10**(i/2.0)))**4)))

constante =  np.mean(I)

fileout = open('s.dat','w')
fileout.write(str(constante))
fileout.close()
