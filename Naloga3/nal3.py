import numpy as np
from matplotlib import pyplot as plt
import scipy
import my_eigensolvers as mygen
from tqdm import tqdm

#np.set_printoptions(edgeitems=30, linewidth=100000, 
#   formatter=dict(float=lambda x: "%.3g" % x))


def expected_q(i:int, j:int):
    """
    Returns <i|q|j> for the regular harmonic oscilator
    """
    if i == j+1: return np.sqrt(i/2)
    if i == j-1: return np.sqrt(j/2)
    else: return 0

def expected_q2(i:int, j:int):
    """
    Returns <i|q^2|j> for the regular harmonic oscilator
    """
    if i == j+2: return np.sqrt((j+1)*(j+2))/2
    if i == j: return (2*j+1)/2
    if i == j-2: return np.sqrt((j)*(j+1))/2
    else: return 0

def expected_q3(i:int, j:int):
    """
    Returns <i|q^3|j> for the regular harmonic oscilator
    """
    if i == j+3: return 1/2**(3/2) * np.sqrt((j+1)*(j+2)*(j+3))
    if i == j+1: return 1/2**(3/2) * 3 * (j+1)**(3/2)
    if i == j-1: return 1/2**(3/2) * 3 * (j)**(3/2)
    if i == j-3: return 1/2**(3/2) * np.sqrt(j*(j-1)*(j-2))
    else: return 0

def expected_q4(i:int, j:int):
    """
    Returns <i|q^4|j> for the regular harmonic oscilator
    """
    if i == j+4: return 1/4 * np.sqrt((j+1)*(j+2)*(j+3)*(j+4))
    if i == j+2: return 1/4 * (3*(j+2)**(3/2)*np.sqrt(j+1) + j * np.sqrt((j+1)*(j+2)))
    if i == j: return 1/4 * 3 * ((j+1)**2 + j**2)
    if i == j-2: return 1/4 * (3*(j-1)**(3/2)*np.sqrt(j) + (j+1) * np.sqrt(j*(j-1)))
    if i == j-4:  return 1/4 * np.sqrt(j*(j-1)*(j-2)*(j-3))
    else: return 0

def expected_q42_iz_navodil(i:int, j:int):
    """
    Returns <i|q^4|j> for the regular harmonic oscilator
    """
    if i == j+4: return 1/2**4 * np.sqrt(2**i * np.math.factorial(i) / 2**j / np.math.factorial(j)) * 1
    if i == j+2: return 1/2**4 * np.sqrt(2**i * np.math.factorial(i) / 2**j / np.math.factorial(j)) * 4 * (2*j+3)
    if i == j: return 1/2**4 * np.sqrt(2**i * np.math.factorial(i) / 2**j / np.math.factorial(j)) * 12 * (2*j**2 +2*j + 1)
    if i == j-2: return 1/2**4 * np.sqrt(2**i * np.math.factorial(i) / 2**j / np.math.factorial(j)) * 16 * j * (2*j**2 - 3*j + 1)
    if i == j-4: return 1/2**4 * np.sqrt(2**i * np.math.factorial(i) / 2**j / np.math.factorial(j)) * 16 * j * (j**3 - 6*j**2 + 11*j - 6)
    else: return 0

"""
#Test, da je vse prav izpeljano
for i in range(1,100):
    for j in range(1,100):
        print(f"{i} {j} ---- {expected_q4(i,j)}, {expected_q42(i,j)}")
"""

def getNondiagonalHamiltonian(N:int, λ:float, q:str) -> np.array:
    """
    N - velikost matrike

    q - parameter, ki pove kateri način izražave <i|q^4|j> uporabimo

    """
    q4 = None
    assert q in ["[q^4]", "[q^2]^2", "[q]^4", "[q^3]^4/3"]
    if q == "[q^4]":
        q4 = expected_q4
    elif q == "[q^2]^2":
        q4 = lambda i, j: expected_q2(i,j)**2
    elif q == "[q]^4":
        q4 = lambda i, j: expected_q(i,j)**4
    elif q == "[q^3]^4/3":
        q4 = lambda i, j: expected_q3(i,j)**(4/3)

    def Hij(i,j):
        return ((i+1/2) + λ*q4(i,j)) if i == j else λ*q4(i,j)

    H = np.zeros((N,N))
    for i in tqdm(range(N), desc = "Building Hamiltonian"):
        for j in range(N):
            H[i,j] = Hij(i,j)

    
    return H


"""
Example
H_nondiagonal = getNondiagonalHamiltonian(100,0.1,"[q^4]")

eigenvalues, eigenvectors = mygen.QR_iteration(H_nondiagonal)

print(eigenvalues)
print(eigenvectors)
"""

def energy_lambda_graph(N:int, q:str, diagonalizationAlgo:str, nlines:int):
    """
    N - velikost matrike

    q - parameter, ki pove kateri način izražave <i|q^4|j> uporabimo

    diagonalizationAlgo {"moj", "np"}

    """
    λs = np.linspace(0,5,100)

    graphs = -np.ones((nlines,len(λs)))
    for i in tqdm(range(len(λs)), desc = "Drawing graph"):
        eigvals, eigvects = None, None
        H = getNondiagonalHamiltonian(N,λs[i],q)
        if diagonalizationAlgo == "moj":
            eigvals, eigvects = mygen.QR_iteration(H)
        else:
            eigvals, eigvects = np.linalg.eigh(H,)

        graphs[:,i] = eigvals[0:nlines]

    for i in range(nlines):
        plt.plot(λs, graphs[i,:])

    plt.xlabel("$\lambda$")
    plt.ylabel("E")
    plt.show()

def energy_N_graph(λ:int, q:str, diagonalizationAlgo:str, nlines:int):
    """
    N - velikost matrike

    q - parameter, ki pove kateri način izražave <i|q^4|j> uporabimo

    diagonalizationAlgo {"moj", "np"}

    """
    Ns = np.arange(4,100,1)



    graphs = -np.ones((nlines,len(Ns)))
    for i in tqdm(range(len(Ns)), desc = "Drawing graph"):
        eigvals, eigvects = None, None
        H = getNondiagonalHamiltonian(Ns[i],λ,q)
        if diagonalizationAlgo == "moj":
            eigvals, eigvects = mygen.QR_iteration(H)
        else:
            eigvals, eigvects = np.linalg.eigh(H)

        graphs[:,i] = eigvals[0:nlines]

    for i in range(nlines):
        plt.plot(Ns, graphs[i,:])

    plt.xlabel("N")
    plt.ylabel("E")
    plt.show()

def energy_N_log_graph(λ:int, q:str, diagonalizationAlgo:str, nlines:int):
    """
    N - velikost matrike

    q - parameter, ki pove kateri način izražave <i|q^4|j> uporabimo

    diagonalizationAlgo {"moj", "np"}

    """
    Ns = np.arange(4,100,1)



    graphs = -np.ones((nlines,len(Ns)))
    for i in tqdm(range(len(Ns)), desc = "Drawing graph"):
        eigvals, eigvects = None, None
        H = getNondiagonalHamiltonian(Ns[i],λ,q)
        if diagonalizationAlgo == "moj":
            eigvals, eigvects = mygen.QR_iteration(H)
        else:
            eigvals, eigvects = np.linalg.eigh(H)

        graphs[:,i] = eigvals[0:nlines]

    for i in range(nlines):
        plt.plot(Ns, np.abs(graphs[i,:]-graphs[i,-1]))

    plt.xlabel("N")
    plt.yscale("log")
    plt.ylabel("E-$E_{\infty}$")
    plt.show()


def plotWaveFunctions(N:int, λ:float, q:str, nlines:int):
    H = getNondiagonalHamiltonian(N,λ,q)
    eigvals, eigvects = np.linalg.eigh(H)

    npoints = 1000
    xs = np.linspace(-5,5, npoints)

    basefuncts = np.zeros((N, npoints))

    haveLoaded = False
    if N == 120:
        try:
            basefuncts = np.load("Naloga3\\basefunctsN=120.npy")
            haveLoaded = True
        except:
            pass

    if haveLoaded == False:
        for i in tqdm(range(N), desc = "Calculating base functions"):
            fi = lambda x: (2**i* np.math.factorial(i) * np.sqrt(np.pi))**(-1/2) * np.exp(-x**2/2) * scipy.special.hermite(i)(x)

            for j in range(npoints):
                basefuncts[i,j] = fi(xs[j])

    #np.save("basefuncts", basefuncts)

    data = np.zeros((nlines, npoints))
    for i in tqdm(range(nlines), desc = "Exbanding eigenstates over base"):
        for j in range(npoints):
            psi = np.zeros(npoints)
            for k in range(N):
                psi += eigvects[k,i] * basefuncts[k,:]
        data[i,:] = psi


    for i in tqdm(range(nlines), desc = "Plotting"):
        plot, = plt.plot(xs, data[i,:]+3*i, zorder = 1)
        plt.axhline(y=3*i, color = plot.get_color(), linewidth = 1)
        plt.fill_between(xs, data[i,:]+3*i, 3*i, color = plot.get_color(), alpha = 0.1)


    plt.plot(xs,(np.square(xs) + λ*np.power(xs, 4)), zorder = 0, linestyle = "dashed", color = "gray")
    plt.ylim(0, 3 * (nlines+1))
    plt.tick_params(
    axis='x',          # changes apply to the x-axis
    which='both',      # both major and minor ticks are affected
    bottom=False,      # ticks along the bottom edge are off
    top=False,         # ticks along the top edge are off
    labelbottom=False) # labels along the bottom edge are off
    plt.tick_params(
    axis='y',          # changes apply to the x-axis
    which='both',      # both major and minor ticks are affected
    left=False,      # ticks along the bottom edge are off
    right=False,         # ticks along the top edge are off
    labelleft=False) # labels along the bottom edge are off
    plt.show()

#N = 120; λ = 1
#plotWaveFunctions(N,λ,"[q]^4", 8)
#plotWaveFunctions(N,λ,"[q^2]^2", 8)
#plotWaveFunctions(N,λ,"[q^3]^4/3", 8)
#plotWaveFunctions(N,λ,"[q^4]", 8)

#energy_lambda_graph(1000,"[q^4]","np", 4)
#energy_lambda_graph(1000,"[q]^4","np", 4)
#energy_lambda_graph(1000,"[q^2]^2","np", 4)
#energy_lambda_graph(1000,"[q^3]^4/3","np", 4)

#energy_N_graph(1,"[q^4]","np", 4)
#energy_N_graph(1,"[q]^4","np", 4)
#energy_N_graph(1,"[q^2]^2","np", 4)
#energy_N_graph(1,"[q^3]^4/3","np", 4)

energy_N_log_graph(1,"[q^4]","np", 4)
energy_N_log_graph(1,"[q^2]^2","np", 4)

    