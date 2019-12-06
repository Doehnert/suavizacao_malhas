from pylab import *

def seidel2d(A,B,Nx,Ny):
    dim = np.size(A,1)
    T = np.zeros(dim)
    n = 100
    for i in range(0,dim):
        T[i] = 0

    for _ in range(0,n):
        for j in range(1,Ny+1):
            for i in range(1,Nx+1):
                aW = 0
                aE = 0
                aS = 0
                aN = 0
                aNE = 0
                aSE = 0
                aNW = 0
                aSW = 0
                P = (i+(j-1)*Nx)-1
                #condicoes de contorno
                aP = A[P][P]
                if j==1 and 1<i<Nx:
                    aN = -A[P][P+Nx]
                    aE = -A[P][P+1]
                    aW = -A[P][P-1]
                    aNE = -A[P][P+Nx+1]
                    aNW = -A[P][P+Nx-1]
                if j==Ny and 1<i<Nx:
                    aS = -A[P][P-(Nx-1)-1]
                    aE = -A[P][P+1]
                    aW = -A[P][P-1]
                    aSE = -A[P][P-(Nx-1)]
                    aSW = -A[P][P-(Nx-1)-2]
                if i==1 and 1<j<Ny:
                    aE = -A[P][P+1]
                    aS = -A[P][P-(Nx-1)-1]
                    aN = -A[P][P+Nx]
                    aNE = -A[P][P+Nx+1]
                    aSE = -A[P][P-(Nx-1)]
                if i==Nx and 1<j<Ny:
                    aW = -A[P][P-1]
                    aS = -A[P][P-(Nx-1)-1]
                    aN = -A[P][P+Nx]
                    aNW = -A[P][P+Nx-1]
                    aSW = -A[P][P-(Nx-1)-2]
                if (i>1 and i<Nx and j>1 and j<Ny):
                    aW = -A[P][P-1]
                    aE = -A[P][P+1]
                    aS = -A[P][P-(Nx-1)-1]
                    aN = -A[P][P+Nx]
                    aNE = -A[P][P+Nx+1]
                    aSE = -A[P][P-(Nx-1)]
                    aNW = -A[P][P+Nx-1]
                    aSW = -A[P][P-(Nx-1)-2]

                if (i<Nx and j<Ny):
                    TNE = T[P+Nx+1]
                else:
                    TNE = 0
                if (i<Nx and j>0):
                    TSE = T[P-(Nx-1)]
                else:
                    TSE = 0
                if (i>0 and j<Ny):
                    TNW = T[P+Nx-1]
                else:
                    TNW = 0
                if (i>0 and j>0):
                    TSW = T[P-(Nx-1)-2]
                else:
                    TSW = 0
                if (i>1):
                    TW = T[P-1]
                else:
                    TW = 0
                if (i<Nx):
                    TE = T[P+1]
                else:
                    TE = 0
                if (j>1):
                    TS = T[P-(Nx-1)-1]
                else:
                    TS = 0
                if (j<Ny):
                    TN = T[P+Nx]
                else:
                    TN = 0

                bP = B[P]
                aP = A[P][P]

                T[P] = ((aW*TW)+(aE*TE)+(aS*TS)+(aN*TN)+(aNE*TNE)+(aSE*TSE)+(aNW*TNW)+(aSW*TSW)+bP)/aP
    return(T)