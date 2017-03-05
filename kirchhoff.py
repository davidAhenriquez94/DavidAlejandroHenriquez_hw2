import numpy as np

datos = []
filein = open('circuito.dat','r')
for i in range(5):
    datos.append(filein.readline())
filein.close()

def crear_matriz(N,R1,R2):
    A = np.zeros((N,N))
    if N == 1:
        A[0,0] = 2.0*R1 + 1.0*R2
        A[0,1] = -1.0*R2
    elif N == 2:
        A[0,0] = 2.0*R1 + 1.0*R2
        A[0,1] = -1.0*R2
        A[1,0] = 1.0*R2
        A[1,1] = 1.0*R1
    else :
        A[0,0] = 2.0*R1 + 1.0*R2
        A[0,1] = -1.0*R2

        for i in range(N-2):
        
            A[i+1,-N+i] = 1.0*R2
            A[i+1,-N+1+i] = 1.0*R1
            A[i+1,-N+2+i] = -1.0*R1
    
            if N%2 == 0:
                A[N-1,N-2] = 1.0*R2
                A[N-1,N-1] = 1.0*R1
            else:
                A[N-1,N-2] = 1.0*R1
                A[N-1,N-1] = 1.0*R2
    
    return A

def resolver_circuito(N,R1,R2,V1,V2):
     A1 = crear_matriz(N,R1,R2)
     A = A1.astype(np.float)
     B1 = np.zeros((N,N))
     B = B1.astype(np.float)
     C1 = np.ones(N)*(1.0*V1+1.0*V2)
     C = C1.astype(np.float)
     D1 = np.zeros(N)
     D = D1.astype(np.float)
     S1 = np.zeros(N)
     S = S1.astype(np.float)

     if N==1:
         S[0] = (V1+V2)/(2.0*R1 + R2)
     elif N==2:
         S[0] = (V1 + V2 + (R2/R1)*(V1+V2))/(2.0*R1 + R2 + (R2**2)/R1)
         S[1] = (V1 + V2 -S[0]*R2)/R1
     
     else:
         B[0,::] = A[0,::]
         B[1,::] = (A[1,0]/A[0,0])*A[0,::]-A[1,::]
         D[1] = (A[1,0]/A[0,0])*C[0]- C[1]

         for i in range(2,N):
             B[i,::] = ((A[i,i-1])/B[i-1,i-1])*B[i-1,::] - A[i,::]
             D[i] = (A[i,i-1]/B[i-1,i-1])*D[i-1] - C[i]
    
         S[-1] = D[-1]*(1.0/B[N-1,N-1])
         for i in range(2,N+1):
             S[-i] = (D[-i] - B[N-i,N-i+1]*S[1-i])*(1.0/B[N-i,N-i])
         
     return S
         
fileout = open('corrientes.dat','w')
sol = resolver_circuito(int(datos[0]),float(datos[1]),float(datos[2]),float(datos[3]),float(datos[4]))
for i in range(sol.size):
    fileout.write("%f\n" %sol[i])
fileout.close()

