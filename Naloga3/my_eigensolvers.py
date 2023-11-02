import numpy as np
from matplotlib import pyplot as plt
from tqdm import tqdm

def QR_Householder(A:np.array) -> tuple:
    A = np.copy(A)

    def Housholder(B:np.array) -> tuple:
        B = np.copy(B)
        N = B.shape[0]

        c1 = B[:,0]

        if np.linalg.norm(c1[1:]) == 0:
            return np.eye(N)

        u = c1 + np.sign(c1[0])* np.linalg.norm(c1) * np.array([1 if i == 0 else 0 for i in range(N)])

        Q = np.eye(N) - 2 * np.outer(u,u) / np.dot(u,u)

        return Q

    Q = np.eye(A.shape[0])
    for i in range(A.shape[0]-1):
        Qi = np.eye(A.shape[0])

        Qi[i:, i:] = Housholder(A[i:,i:])

        Q = Q.dot(Qi)
        A = Qi.dot(A)

    return (A, Q)

def Tridiag_Householder(A:np.array) -> tuple:
    A = np.copy(A)

    def Housholder(B:np.array) -> np.array:
        B = np.copy(B)
        N = B.shape[0]

        c1 = B[:,0]

        if np.linalg.norm(c1[1:]) == 0:
            return np.eye(N)

        u = c1 + np.sign(c1[0])* np.linalg.norm(c1) * np.array([1 if i == 0 else 0 for i in range(N)])

        Q = np.eye(N) - 2 * np.outer(u,u) / np.dot(u,u)

        return Q

    Q = np.eye(A.shape[0])
    for i in tqdm(range(1, A.shape[0]-1), desc="Tridiag"):
        Qi = np.eye(A.shape[0])
        Qi[i:, i:] = Housholder(A[i:,i-1:])

        Q = Q.dot(Qi)
        A = Qi.dot(A).dot(Qi)

    return (A, Q)

def QR_iteration(M:np.array, do_tridiag:bool = True) -> tuple:
    if do_tridiag:
        T, Q = Tridiag_Householder(M)
    else:
        T = np.copy(M)
        Q = np.eye(T.shape[0])

    def offDiag(A:np.array) -> float:
        return np.sum(A) - np.trace(A)
    
    epsilon = 10**-10

    A = np.copy(T)
    Z = np.eye(A.shape[0])
    while abs(offDiag(A)) > epsilon:
        print(f"{epsilon} {offDiag(A)}", end="\r")
        Ri,Qi = QR_Householder(A)
        A = Ri.dot(Qi)
        Z = Z.dot(Qi)

    eigvects = Q.dot(Z).dot(Q.T)
    eigvals = np.diag(A)

    inds = np.argsort(eigvals)
    

    return (eigvals[inds], eigvects[:,inds])

import timeit

def QR_iteration_comparison():
    start = timeit.default_timer()
    eigvals, eigvects = QR_iteration(np.array([[1,1,1], [1,0,1], [1,1,2]]))
    
    print(f"Time taken is {timeit.default_timer()-start}")
    print(eigvals)
    print(eigvects)
    print("\n")

    start = timeit.default_timer()
    eigvals, eigvects = np.linalg.eigh(np.array([[1,1,1], [1,0,1], [1,1,2]]))

    print(f"Time taken is {timeit.default_timer()-start}")
    print(eigvals)
    print(eigvects)

def getRandomMatrix(N:int):
    return np.random.randint(-10,10,size=(N,N))

def getRandomSymetricMatrix(N:int):
    b = np.random.randint(-10,10,size=(N,N))
    return (b + b.T)/2


Ns = [1,2,3,4,5,6,7,8,9,10,15,25,30]

y1,y2,y3,y4 = [],[],[],[]
for N in Ns:
    M = getRandomSymetricMatrix(N)
    print(N)
    for i in range(3):
        y = []
        start = timeit.default_timer() * 1000
        eigvals, eigvects = QR_iteration(M, do_tridiag=True)
        print(f"Time taken is {timeit.default_timer()*1000-start}")
        y.append(timeit.default_timer()*1000-start)
    y1.append(np.average(y))

    for i in range(3):
        y = []
        start = timeit.default_timer() * 1000
        eigvals, eigvects = QR_iteration(M, do_tridiag=False)
        print(f"Time taken is {timeit.default_timer()*1000-start}")
        y.append(timeit.default_timer()*1000-start)
    y2.append(np.average(y))

Ns2 = np.arange(1,1000)
for N in Ns2:
    M = getRandomSymetricMatrix(N)
    print(N)

    for i in range(10):
        y = []
        start = timeit.default_timer() * 1000
        eigvals, eigvects = np.linalg.eigh(M)
        print(f"{N} Time taken is {timeit.default_timer()*1000-start}")
        y.append(timeit.default_timer()*1000-start)
    y3.append(np.average(y))

    for i in range(10):
        y = []
        start = timeit.default_timer() * 1000
        eigvals, eigvects = np.linalg.eig(M)
        print(f"{N} Time taken is {timeit.default_timer()*1000-start}")
        y.append(timeit.default_timer()*1000-start)
    y4.append( np.average(y))

plt.plot(Ns,y1, label="QR z tridiaggonalizacijo")
plt.plot(Ns,y2, label="QR")
plt.plot(Ns2,y3, label="numpy.eigh")
plt.plot(Ns2,y4, label="numpy.eig")
plt.ylabel("Čas izvajanja [ms]")
plt.xlabel("N")
plt.legend()
plt.show()


plt.plot(Ns,np.array(y1)/1000, label="QR z tridiaggonalizacijo")
plt.plot(Ns,y2/1000, label="QR")
plt.ylabel("Čas izvajanja [s]")
plt.xlabel("N")
plt.legend()
plt.show()


plt.plot(Ns2,np.array(y1)/y3, label="numpy.eigh")
plt.plot(Ns2,y4, label="numpy.eig")
plt.ylabel("Čas izvajanja [ms]")
plt.xlabel("N")
plt.legend()
plt.show()





"""
M = getRandomSymetricMatrix(2)

start = timeit.default_timer() * 1000
eigvals, eigvects = QR_iteration(M, do_tridiag=True)
print(f"Time taken is {timeit.default_timer()*1000-start}")

start = timeit.default_timer() * 1000
eigvals, eigvects = QR_iteration(M, do_tridiag=False)
print(f"Time taken is {timeit.default_timer()*1000-start}")

start = timeit.default_timer()*1000
eigvals, eigvects = np.linalg.eigh(M)
print(f"Time taken is {timeit.default_timer()*1000-start}")

start = timeit.default_timer()*1000
eigvals, eigvects = np.linalg.eig(M)
print(f"Time taken is {timeit.default_timer()*1000-start}")
"""