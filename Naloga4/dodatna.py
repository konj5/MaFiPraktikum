import numpy as np
import matplotlib.pyplot as plt
import scipy

def scientificToFloat(s: str):
    a, b = s.split("E")
    a = float(a); b = float(b)

    return a * 10**b

for freq in [300, 600, 900, 933]:
    with open(f"Naloga4\\AkRes\\Zvočni profil {freq}Hz.txt", "r") as f:
        lines = []

        for line in f.readlines():
            line = line.strip()
            lines.append(line.split("\t"))


        amp = []
        dev = []

        for line in lines:
            amp.append(scientificToFloat(line[0]))
            dev.append(scientificToFloat(line[1]))


        odmik = [x for x in range(0,len(amp),1)]

    dx = odmik[1]-odmik[0]
    N = len(odmik)
    sampling_frequency = 1/dx
    nuc = sampling_frequency/2
    nu = np.linspace(-nuc,nuc, N, endpoint=False)

    fourier = np.fft.fft(dev)/N
    fourier = np.roll(fourier, N//2)

    fig, axs = plt.subplots(3,1)

    axs[0].plot(odmik, dev)
    axs[0].scatter(odmik, dev, color = "black", s = 6)
    axs[0].set_xlabel("odmik")
    axs[0].set_ylabel("Amplituda")

    axs[1].plot(nu, np.real(fourier))
    axs[1].scatter(nu,  np.real(fourier), color = "black", s = 6)
    axs[1].set_xlabel("k")
    axs[1].set_ylabel("Re")

    axs[2].plot(nu, np.imag(fourier))
    axs[2].scatter(nu, np.imag(fourier), color = "black", s = 6)
    axs[2].set_xlabel("k")
    axs[2].set_ylabel("Im")

    fig.tight_layout()

    plt.show()

