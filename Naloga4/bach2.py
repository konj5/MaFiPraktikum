import numpy as np
from matplotlib import pyplot as plt
import scipy

use_fft = False


fnames = ["Naloga4\\Bach.44100.txt","Naloga4\\Bach.11025.txt","Naloga4\\Bach.5512.txt","Naloga4\\Bach.2756.txt","Naloga4\\Bach.1378.txt","Naloga4\\Bach.882.txt"]
fnames.reverse()

if use_fft == True:
    fnames = ["Naloga4\\data " + fname[13:-4] + " FFT" for fname in fnames]
else:
    fnames = ["Naloga4\\data " + fname[13:-4] + " DFT" for fname in fnames]

i = 0
fig, axs = plt.subplots(6,1)

for fname in fnames:
    vk, tn, signal, fourier = [],[],[],[]
    with open(fname, "r") as f:
        lines = f.readlines()

        for line in lines:
            line = line[1:-1]
            split = line.split(",")
            vk.append(float(split[0]))
            tn.append(float(split[1]))
            signal.append(float(split[2]))
            fourier.append(complex(split[3][1:-1]))

    if i <= 4:    
        axs[i].plot(vk, np.abs(fourier))
        axs[i].set_xlabel("$\\nu$")
        axs[i].set_ylabel("$|F($\\nu$)|$")
        axs[i].set_title(f"{fname[13:-4]}Hz")
        axs[i].set_xlim((-10000, 10000))
        i += 1
    
    elif i == 5:
        axs[i].plot(vk, np.abs(fourier))
        axs[i].set_xlabel("$\\nu$")
        axs[i].set_ylabel("$|F($\\nu$)|$")
        axs[i].set_title(f"{fname[13:-4]}Hz")
        axs[i].set_xlim((-10000, 10000))
        i += 1

        for j in range(6):
            axs[j].fill(vk, np.abs(fourier), color = "orange")


plt.show()








