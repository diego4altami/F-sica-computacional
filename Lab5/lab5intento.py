import numpy as np
import numpy.random as nr

def permutation():
    global order,N
    order = nr.permutation(N)
def vecinocuad(s,j):
    global L,empty

    if j == 0: #izquierda
        if s%L == L-1:
            return empty
        else:
            return s+1
        
    elif j==1: #arriba
        if s//L == L-1:
            return empty
        else:
            return s+L
        
    elif j==2: #derecha
        if s%L == 0:
            return empty
        else:
            return s-1
        
    elif j==3: #abajo
        if s//L == 0:
            return empty
        else:
            return s-L
def boundaries():
    global N,nn
    for s in range(N):
        for j in range(4):
            nn[s,j] = vecinocuad(s,j)
def findroot(r):
    global parent
    while parent[r] >= 0:
        r = parent[r]
    return r
def mergeroots(r1,r2):
    global spanclussize,parent,touchesLeft,touchesRight
    if r1 == r2:
        return r1
    elif -parent[r1] > -parent[r2]:
        parent[r1] += parent[r2]
        parent[r2] = r1
        touchesLeft[r1] = touchesLeft[r1] or touchesLeft[r2] 
        touchesRight[r1] = touchesRight[r1] or touchesRight[r2]
        if touchesLeft[r1] and touchesRight[r1]:
            spanclussize = -parent[r1]
        return r1 
    else:
        parent[r2] += parent[r1]
        parent[r1] = r2
        touchesLeft[r2] = touchesLeft[r2] or touchesLeft[r1]
        touchesRight[r2] = touchesRight[r2] or touchesRight[r1]
        if touchesLeft[r2] and touchesRight[r2]:
            spanclussize = -parent[r2]
        return r2  
def percolate():
    global spanclussize,order,parent,nn
    BIG = []
    PSpan = []
    big=0
    spanclussize = 0
    
    for i in range(N): parent[i]=empty    

    for i in range(N):
        r1=s1=order[i]
        parent[s1] = -1
                
        for j in range(4):
            s2=nn[s1, j]
            if s2 != empty:
                if parent[s2] != empty:
                    r1 = mergeroots(r1,findroot(s2))
                    if -parent[r1]>big: big=-parent[r1]
        BIG.append(big)
        PSpan.append(spanclussize/(i+1.))

    return np.array(BIG), np.array(PSpan)
global N, empty, parent, nn, order, spanclussize,touchesLeft,touchesRight

L = 128
N = L**2
empty=-(N+1)
muestras=100

def multiplepercolate(L):
    global N,muestras
    bigprom=np.zeros(N)
    Pspanprom=np.zeros(N)
    Probprom=np.zeros(N)
    nn=np.zeros((N, 4), dtype=int)
    boundaries() 

    for i in range (muestras):
        order=np.zeros(N, dtype=int) 
        parent=np.zeros(N, dtype=int) 
        touchesLeft =  [i%L == 0 for i in range(N)]
        touchesRight = [i%L == L-1 for i in range(N)]

        permutation()    
        BIG,PSpan=percolate()
        Prob = np.copy(PSpan)
        Prob[PSpan>0]=1

        bigprom+=BIG
        Pspanprom+=PSpan
        Probprom+=Prob

    bigprom/=muestras
    Pspanprom/=muestras 
    Probprom/=muestras

    return bigprom,Pspanprom,Probprom

bigprom,Pspanprom,Probprom=multiplepercolate()

print("El arreglo más grande promedio como funció de p es: ",bigprom)
print("El Pspan promedio como funció de p es: ",Pspanprom)
print("La probabilidad de que percole promedio como funció de p es: ",Probprom)

size=np.array([16,32,64,128,256])