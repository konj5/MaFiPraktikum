import numpy as np
from matplotlib import pyplot as plt
import scipy
from tqdm import tqdm

def simple_fourier_transform(fn):
    fn = np.array(fn)
    N = len(fn)
    Fk = np.zeros(N, dtype=complex)

    for k in tqdm(range(N), desc="fourier trasform"):
        sum = 0
        for n in range(N):    
            sum += fn[n] * np.exp(-2j * np.pi * k * n / N)
        Fk[k] = sum
    return Fk


def simple_inverse_fourier_transform(Fk):
    Fk = np.array(Fk)
    N = len(Fk)
    fn = np.zeros(N, dtype=complex)

    for n in range(N):
        sum = 0
        for k in range(N):    
            sum += Fk[k] * np.exp(2j * np.pi * k * n / N)
        fn[n] = sum
    return fn

def getFrequencySpace(fn,T):
    N = len(fn)
    dt = T/N
    nuc = 0.5/dt
    vk= np.linspace(-nuc,nuc,N,endpoint=False)
    return vk


""" Od tu naprej so samo na stavitve za risanje raznih grafov"""

"""
N = 200
T = 100
aliasing = False
aliasingnu = None
modulation = False
sigma = 7

#f = lambda t: 1
#f = lambda t: np.sin(2*np.pi * 20/T *t) + np.cos(2*np.pi* 2/T *t)
f = lambda t: np.cos(2*np.pi* 120/T *t); aliasing = True; aliasingnu = 0.8
#f = lambda t: 1/np.sqrt(2*np.pi*sigma**2) * np.exp(-1/2 / sigma**2 * (t-T/2)**2); modulation = True
#f = lambda t: 1/np.sqrt(2*np.pi*sigma**2) * np.exp(-1/2 / sigma**2 * (t)**2) if t < T/2 else 1/np.sqrt(2*np.pi*sigma**2) * np.exp(-1/2 / sigma**2 * (T-t)**2)
#f = lambda t: t/T*2 if t < T/2 else 2 - t/T*2
#f = lambda t: np.exp(-t**2/50)

tn = np.linspace(0,T,N,endpoint=False)
#tn = np.linspace(0,T,N,endpoint=True) #Če želimo videti, kako hudo se pokvari periodičnost
fn = [f(t) for t in tn]
Fk = simple_fourier_transform(fn) / N
Fk_rolled = np.roll(Fk, int(N/2))

fn2 = simple_inverse_fourier_transform(Fk)


dt = T/N
nuc = 0.5/dt
vk_unrolled = np.linspace(0,2*nuc,N,endpoint=False)
vk= np.linspace(-nuc,nuc,N,endpoint=False)

"""

"""
if aliasing == False and modulation == False:
    fig, axs = plt.subplots(3,1)
    axs[0].plot(np.linspace(0,T,10*N), [f(t) for t in np.linspace(0,T,10*N)], zorder = 0)
    axs[0].scatter(tn, fn, color = "black", s = 6, zorder = 1)
    axs[0].set_xlabel("t")
    axs[0].set_ylabel("f(t)")

    axs[1].plot(vk, np.real(Fk_rolled))
    axs[1].set_xlabel("$\\nu$")
    axs[1].set_ylabel("Re[F($\\nu$)]")

    axs[2].plot(vk, np.imag(Fk_rolled))
    axs[2].set_xlabel("$\\nu$")
    axs[2].set_ylabel("Im[F($\\nu$)]")
    fig.tight_layout()
    plt.show()

if aliasing == True:
    fig, axs = plt.subplots(3,1)
    axs[0].plot(np.linspace(0,T,10*N)[0:250], [f(t) for t in np.linspace(0,T,10*N)][0:250], zorder = 0)

    axs[0].plot(np.linspace(0,T,10*N)[0:250], [1/2 * np.cos(2*np.pi* aliasingnu *t) + 1/2*np.cos(2*np.pi* -aliasingnu *t) for t in np.linspace(0,T,10*N)][0:250], zorder = 1, color = "red", linestyle = "dashed")

    axs[0].scatter(tn[0:25], fn[0:25], color = "black", s = 6, zorder = 2)
    axs[0].set_xlabel("t")
    axs[0].set_ylabel("f(t)")

    axs[1].plot(vk, np.real(Fk_rolled))
    axs[1].set_xlabel("$\\nu$")
    axs[1].set_ylabel("Re[F($\\nu$)]")

    axs[2].plot(vk, np.imag(Fk_rolled))
    axs[2].set_xlabel("$\\nu$")
    axs[2].set_ylabel("Im[F($\\nu$)]")
    fig.tight_layout()
    plt.show()

if modulation == True:
    fig, axs = plt.subplots(3,1)
    axs[0].plot(np.linspace(0,T,10*N), [f(t) for t in np.linspace(0,T,10*N)], zorder = 0)
    axs[0].scatter(tn, fn, color = "black", s = 6, zorder = 1)
    axs[0].set_xlabel("t")
    axs[0].set_ylabel("f(t)")

    axs[1].plot(vk, np.real(Fk_rolled), label = "brez popravka")

    Fk_mod = np.array([Fk_rolled[i]  *np.exp(2j* np.pi * vk[i] * T/2) for i in range(len(Fk))])
    axs[1].plot(vk, np.real(Fk_mod), label = "z popravkom")
    axs[1].set_xlabel("$\\nu$")
    axs[1].set_ylabel("Re[F($\\nu$)]")
    axs[1].legend()

    axs[2].plot(vk, np.imag(Fk_rolled))
    axs[2].set_xlabel("$\\nu$")
    axs[2].set_ylabel("Im[F($\\nu$)]")
    fig.tight_layout()
    plt.show()

"""
"""
fig, axs = plt.subplots(5,1)
axs[0].plot(np.linspace(0,T,10*N), [f(t) for t in np.linspace(0,T,10*N)], zorder = 0)
axs[0].scatter(tn, fn, color = "black", s = 6, zorder = 1)
axs[0].set_xlabel("t")
axs[0].set_ylabel("f(t)")

axs[1].plot(vk, np.real(Fk_rolled))
axs[1].set_xlabel("$\\nu$")
axs[1].set_ylabel("Re[F($\\nu$)]")

axs[2].plot(vk, np.real(fn2))
axs[2].set_xlabel("t")
axs[2].set_ylabel("Re[$f_{inverz}$(t)]")

axs[3].plot(vk, np.imag(fn2))
axs[3].set_xlabel("t")
axs[3].set_ylabel("im[$f_{inverz}$(t)]")

axs[4].plot(vk, np.abs(fn2-fn) )
axs[4].set_xlabel("t")
axs[4].set_ylabel("Absoutna napaka")
axs[4].set_yscale("log")

fig.tight_layout()
plt.show()
"""
"""
import timeit

Ns = np.arange(1,1000)
times = np.zeros(Ns.shape[0])
fftimes = np.zeros(Ns.shape[0])

for i in tqdm(range(len(Ns)-1, -1, -1)):
    tn = np.linspace(0,T,Ns[i],endpoint=False)
    fn = [f(t) for t in tn]

    start = timeit.default_timer()
    simple_fourier_transform(fn)
    times[i]=(timeit.default_timer()-start)

    start = timeit.default_timer()
    np.fft.fft(fn)
    fftimes[i]=(timeit.default_timer()-start)

plt.scatter(Ns,times, label = "DFT", s = 5)
plt.scatter(Ns,fftimes, label = "FFT", s = 5)
plt.legend()
plt.xlabel("N")
plt.ylabel("Čas izvajanja")
plt.yscale("log")
plt.show()
"""

"""
xs = np.linspace(0,2,200, endpoint=False)
ys = []

skiped = 0
i = 0
while i < len(xs):
    if xs[i] == 1:
        skiped += 32
        print("jit2")

    
    ys.append(np.sin(2*np.pi* 6 *(xs[i] + skiped * 2/200)) + np.cos(2*np.pi*(xs[i] + skiped * 2/200)))
    i+=1

plt.plot(xs,ys)
plt.show()
"""
"""
xs = np.linspace(0,400,200, endpoint=False)
sigma = 20
ys = [1/np.sqrt(2*np.pi*sigma**2) * np.exp(-1/2 / sigma**2 * (np.abs(t) % 100)**2) for t in xs]

plt.plot(xs,ys)
plt.show()
"""
"""
xs = np.linspace(0,800,200, endpoint=False)
sigma = 20
ys = [1/np.sqrt(2*np.pi*sigma**2) * np.exp(-1/2 / sigma**2 * ((np.abs(t) % 200)-100)**2) for t in xs]

plt.plot(xs,ys)
plt.show()
"""