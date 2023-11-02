import numpy as np
from matplotlib import pyplot as plt
import scipy
from tqdm import tqdm
from MyDFT import simple_fourier_transform as DFT


use_fft = False


fnames = ["Naloga4\\Bach.44100.txt","Naloga4\\Bach.11025.txt","Naloga4\\Bach.5512.txt","Naloga4\\Bach.2756.txt","Naloga4\\Bach.1378.txt","Naloga4\\Bach.882.txt"]
fnames.reverse()

vks, tns, sigs, fours = [],[],[],[]

readfiles = True
if readfiles == False:
    for i in tqdm(range(len(fnames)), desc="overall progress"):

        with open(fnames[i], mode="r") as f:
            lines = f.readlines()

            signal = []
            for line in lines:
                signal.append(float(line))

            signal = np.asarray(signal)



        N = len(signal)
        sampling_frequency = float(fnames[i][13:-4])
        print(sampling_frequency)


        nuc = sampling_frequency / 2

        if use_fft == False:
            fourier = DFT(signal)/N
        else:
            fourier = np.fft.fft(signal)/N

        fourier = np.roll(fourier, len(fourier)//2)


        vk = np.linspace(-nuc,nuc,N,endpoint=False)
        tn = np.array([i * 1/sampling_frequency for i in range(N)])

        if use_fft == False:
            with open(f"Naloga4\\data {fnames[i][13:-4]} DFT", mode="w") as f:
                for i in range(len(vk)):
                    f.write(f"{vk[i],tn[i],signal[i],fourier[i]}\n")
        else:
            with open(f"Naloga4\\data {fnames[i][13:-4]} FFT", mode="w") as f:
                for i in range(len(vk)):
                    f.write(f"{vk[i],tn[i],signal[i],fourier[i]}\n")

        vks.append(vk)
        tns.append(tn)
        sigs.append(signal)
        fours.append(fourier) 


    for i in range(len(fnames)):

        fig, axs = plt.subplots(3,1)
        axs[0].plot(tns[i], sigs[i])
        axs[0].set_xlabel("t")
        axs[0].set_ylabel("f(t)")

        axs[1].plot(vks[i], np.real(fours[i]))
        axs[1].set_xlabel("$\\nu$")
        axs[1].set_ylabel("Re[F($\\nu$)]")

        axs[2].plot(vks[i], np.imag(fours[i]))
        axs[2].set_xlabel("$\\nu$")
        axs[2].set_ylabel("Im[F($\\nu$)]")
        fig.suptitle(f"Fourierova analiza {fnames[i][13:]}")
        fig.tight_layout()
        plt.show()

else:
    if use_fft == True:
        fnames = ["Naloga4\\data " + fname[13:-4] + " FFT" for fname in fnames]
    else:
        fnames = ["Naloga4\\data " + fname[13:-4] + " DFT" for fname in fnames]

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


        fig, axs = plt.subplots(2,1)

        axs[0].plot(vk, np.real(fourier))
        axs[0].set_xlabel("$\\nu$")
        axs[0].set_ylabel("Re[F($\\nu$)]")

        axs[1].plot(vk, np.imag(fourier))
        axs[1].set_xlabel("$\\nu$")
        axs[1].set_ylabel("Im[F($\\nu$)]")

        fig.suptitle(f"Fourierova analiza vzorca vzorčenega pri {fname[13:-4]} Hz")
        fig.tight_layout()
        plt.show()