import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as ani

def non_FT_Corelation(f:list,g:list,n:int):
    N = len(g)
    sum = 0
    for i in range(N):
        try:
            sum += g[i+n] * f[i]
        except IndexError:
            break

    return sum / (N-n)

def rescaleAutocorelation(autocorr, h):
    avg = np.average(h)
    autocorr = np.asarray(autocorr)

    return (autocorr - avg**2) / (autocorr[0]-avg**2)

def corelation_FFT(f:list,g:list):
    if len(f) > len(g):
        g = np.concatenate((g, np.zeros(len(f)-len(g))))
    elif len(f) < len(g):
        f = np.concatenate((f, np.zeros(len(g)-len(f))))



    f = np.asarray(f); g = np.asarray(g)
    N = len(f)
    f = np.concatenate((f, np.zeros(N)))
    g = np.concatenate((g, np.zeros(N)))

    ns = np.arange(-N,N)

    F = np.fft.fft(f)/N
    G = np.fft.fft(g)/N

    temp = np.conjugate(F) * G
    temp = np.fft.ifft(temp)

    return ns, np.real(np.fft.ifftshift(temp))

def BADcorelation_FFT_DO_NOT_USE_UNLESS_USED_AS_A_BAD_EXAMPLE(f,g):
    assert len(f) == len(g)
    f = np.asarray(f); g = np.asarray(g)
    N = len(f)

    ns = np.arange(-N//2,N//2)

    F = np.fft.fft(f)/N
    G = np.fft.fft(g)/N

    temp = np.conjugate(F) * G
    temp = np.fft.ifft(temp)

    return ns, np.real(np.fft.ifftshift(temp))

def BoxFunction(x,x0,x1):
    if x >= x0 and x <= x1:
        return 1
    
    return 0


def AnimateCorelation(f,g):
    assert len(f) == len(g)
    speed = 3

    fig, axs = plt.subplots(2,1)
    
    xs = np.arange(0,len(f))
    tau = np.arange(-len(f), len(f), speed, dtype=int)
    ns, corr = corelation_FFT(f,g)

    corr = corr / np.max(corr)


    fplot = axs[0].plot(xs,f,label=f"f(x-$\\tau$)")
    axs[0].plot(xs,g,label= "g(x)")
    axs[0].legend()

    axs[1].plot(ns,corr)
    point = axs[1].scatter(tau[0], corr[0])
    axs[1].set_title("Korelacija (Normirana)")
    axs[1].set_xlabel("$\\tau$")


    def animate(frame, fplot):
        f_t = []
        for i in range(len(f)):
            try:
                if i-tau[frame] >= 0:
                    f_t.append(f[i-tau[frame]])
                else:
                    f_t.append(0)
            except IndexError:
                f_t.append(0)

        
        fplot.set_ydata(f_t)
        fplot.set_label("f(x-$\\tau$), $\\tau = {tau[frame]}$")
        point.set_offsets((tau[frame], corr[speed*frame]))
        return (fplot, point)


    anim = ani.FuncAnimation(fig=fig, func=animate, fargs=(fplot), frames=2 * len(f)//speed, interval=1)
    fig.tight_layout()

    
    writer = ani.PillowWriter(fps=15,
                                    metadata=dict(artist='Me'),
                                    bitrate=1800)
    
    #anim.save('Corelation.gif', writer=writer)

    plt.show()


def AnimateBADCorelation_DO_NOT_USE_UNLESS_USED_AS_A_BAD_EXAMPLE(f,g):
    assert len(f) == len(g)
    speed = 3

    fig, axs = plt.subplots(2,1)
    
    xs = np.arange(0,len(f))
    tau = np.arange(-len(f)//2, len(f)//2, speed, dtype=int)
    ns, corr = BADcorelation_FFT_DO_NOT_USE_UNLESS_USED_AS_A_BAD_EXAMPLE(f,g)

    corr = corr / np.max(corr)


    fplot = axs[0].plot(xs,f,label=f"f(x-$\\tau$)")
    axs[0].plot(xs,g,label= "g(x)")
    axs[0].legend()

    axs[1].plot(ns,corr)
    point = axs[1].scatter(tau[0], corr[0])
    axs[1].set_title("Korelacija (Normirana)")
    axs[1].set_xlabel("$\\tau$")



    def animate(frame, fplot):
        f_t = np.roll(f,tau[frame])

        
        fplot.set_ydata(f_t)
        fplot.set_label("f(x-$\\tau$), $\\tau = {tau[frame]}$")
        point.set_offsets((tau[frame], corr[speed*frame]))
        return (fplot, point)


    anim = ani.FuncAnimation(fig=fig, func=animate, fargs=(fplot), frames= len(f)//speed, interval=1)
    fig.tight_layout()

    
    writer = ani.PillowWriter(fps=15,
                                    metadata=dict(artist='Me'),
                                    bitrate=1800)
    
    anim.save('BadCorelation.gif', writer=writer)

    plt.show()

"""
f = [np.sin(2*np.pi*x)*BoxFunction(x,0,1) + BoxFunction(x,3,5) + BoxFunction(x,9,10) for x in np.linspace(0,10,1000)]
g = [BoxFunction(x,6,8) for x in np.linspace(0,10,1000)]


fig,axs = plt.subplots(1,2)
axs[0].plot(f,label= "f(x)")
axs[1].plot(g,label= "g(x)")
plt.show()

AnimateCorelation(f,g)
AnimateBADCorelation_DO_NOT_USE_UNLESS_USED_AS_A_BAD_EXAMPLE(f,g)
"""


"""
f = [np.sin(3*np.pi*(x-1)) * BoxFunction(x,1,2) for x in np.linspace(0,5,500)]
g = [BoxFunction(x,0,1) + BoxFunction(x,1.1,1.2) + np.log(x+1)*BoxFunction(x,3,4.6) for x in np.linspace(0,5,500)]

f_with_noise = [x + (np.random.rand()-1/2)*10 for x in f]


fig, axs = plt.subplots(4,2)
fig.set_size_inches(w=11,h=10*0.8)

original1 = f
original2 = g
signal = f_with_noise

axs[0][0].plot(np.linspace(0,2,len(original1)), original1)
axs[0][0].set_xlabel("x")
axs[0][0].set_title("f(x)")

axs[0][1].plot(np.linspace(0,2,len(original2)), original2)
axs[0][1].set_xlabel("x")
axs[0][1].set_title("g(x)")

sign, sigacor = corelation_FFT(signal,signal)

axs[1][0].plot(np.linspace(0,5,len(signal)), signal)
axs[1][0].set_xlabel("x")
axs[1][0].set_title("f(x), z naključnim šumom")

axs[1][1].plot(sign , sigacor)
axs[1][1].set_xlabel("Zamik")
axs[1][1].set_title("Avtokorelacija zašumljenega f(x)")

ans1, acor1 = corelation_FFT(original1, original1)
ans2, acor2 = corelation_FFT(original2, original2)

axs[2][0].plot(ans1, acor1)
axs[2][0].set_xlabel("Zamik")
axs[2][0].set_title("Avtokorelacija f(x)")

axs[2][1].plot(ans2 , acor2)
axs[2][1].set_xlabel("Zamik")
axs[2][1].set_title("Avtokorelacija g(x)")

ns1, corr1 = corelation_FFT(original1, signal)
ns2, corr2 = corelation_FFT(original2, signal)

axs[3][0].plot(ns1, corr1)
axs[3][0].set_xlabel("Zamik")
axs[3][0].set_title("Korelacija f(x) z zašumljenim f(x)")

axs[3][1].plot(ns2, corr2)
axs[3][1].set_xlabel("Zamik")
axs[3][1].set_title("Korelacija g(x) z zašumljenim f(x)")
        

fig.tight_layout()
plt.show()
"""


def getDataFromFile(fname:str):
    with open(fname, mode="r") as f:
        lines = f.readlines()

    data = []
    for line in lines:
        data.append(float(line))
    return data

"""
fnames = ["Naloga5\\bubomono.txt", "Naloga5\\bubo2mono.txt", "Naloga5\\mix.txt", "Naloga5\\mix1.txt", "Naloga5\\mix2.txt", "Naloga5\\mix3.txt"]

for mix_fname in fnames[2:]:

    fig, axs = plt.subplots(4,2)
    fig.set_size_inches(w=11,h=10*0.8)

    original1 = getDataFromFile(fnames[0])[1:]
    original2 = getDataFromFile(fnames[1])[1:]
    signal = getDataFromFile(mix_fname)[1:]

    axs[0][0].plot(np.linspace(0,5,len(original1)), original1)
    axs[0][0].set_xlabel("Čas[s]")
    axs[0][0].set_title("Prva sova")

    axs[0][1].plot(np.linspace(0,5,len(original2)), original2)
    axs[0][1].set_xlabel("Čas[s]")
    axs[0][1].set_title("Druga sova")

    signame = mix_fname.split("\\")[1].split(".")[0]
    sign, sigacor = corelation_FFT(signal,signal)

    axs[1][0].plot(np.linspace(0,5,len(signal)), signal)
    axs[1][0].set_xlabel("Čas[s]")
    axs[1][0].set_title(f"Posnetek {signame}.waw")

    axs[1][1].plot(sign* 5/len(signal), sigacor)
    axs[1][1].set_xlabel("Zamik[s]")
    axs[1][1].set_title(f"Avtokorelacija {signame}.waw")

    ans1, acor1 = corelation_FFT(original1, original1)
    ans2, acor2 = corelation_FFT(original2, original2)

    axs[2][0].plot(ans1 * 5/len(original1) , acor1)
    axs[2][0].set_xlabel("Zamik[s]")
    axs[2][0].set_title("Avtokorelacija prve sove")

    axs[2][1].plot(ans2 * 5/len(original2) , acor2)
    axs[2][1].set_xlabel("Zamik[s]")
    axs[2][1].set_title("Avtokorelacija druge sove")

    ns1, corr1 = corelation_FFT(original1, signal)
    ns2, corr2 = corelation_FFT(original2, signal)

    axs[3][0].plot(ns1 * 5/len(original1) , corr1)
    axs[3][0].set_xlabel("Zamik[s]")
    axs[3][0].set_title(f"Korelacija prve sove z {signame}.waw")

    axs[3][1].plot(ns2 * 5/len(original2) , corr2)
    axs[3][1].set_xlabel("Zamik[s]")
    axs[3][1].set_title(f"Korelacija druge sove z {signame}.waw")
            

    fig.tight_layout()
    plt.show()

"""

