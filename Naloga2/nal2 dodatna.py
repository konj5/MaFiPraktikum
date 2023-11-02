import numpy as np
from matplotlib import pyplot as plt
import random
import statsmodels.api as sm
#from progress.bar import Bar
from tqdm import tqdm



random.seed(69420)

#randomly generates angle from even distribution
def randomPhi(phimin = 0, phimax = np.pi * 2):
    return random.random() * (phimax-phimin) + phimin

#randomly generates step length from pareto distribution
def randomL(Lmin, mu):
    rho = random.random()
    return Lmin * (1-rho)**(1/(1-mu))

#randomly generates step length from pareto distribution
def randomT(Tmin, nu):
    rho = random.random()
    return Tmin * (1-rho)**(1/(1-nu))

#Test for randomL function
"""
Lmin = 1
mu = 4

ls = [randomL(Lmin,mu) for _ in range(10000)]
xs = np.linspace(Lmin, 10,1000)
ys = [(mu-1) * Lmin**(mu-1) * x**(-mu) for x in xs]


plt.hist(ls,bins=1000, density=True)
plt.plot(xs,ys)
plt.show()
"""

def generateLevy(Nsteps, Lmin, mu):
    xpos, ypos, Ls = [0], [0], []
    x, y = 0, 0

    for _ in range(Nsteps):
        #Get and store L and phi
        phi = randomPhi()
        L = randomL(Lmin=Lmin, mu=mu)
        Ls.append(L)

        #Update position        
        x = x + L * np.cos(phi)
        y = y + L * np.sin(phi)
        xpos.append(x)
        ypos.append(y)

    return (xpos,ypos,Ls)

def getLevies(Nruns, Nsteps, Lmin, mu):
    xcoords = np.zeros((Nruns,Nsteps+1))
    ycoords = np.zeros((Nruns,Nsteps+1))
    Lss = np.zeros((Nruns,Nsteps))

    for i in range(Nruns):
        xpos, ypos, Ls = generateLevy(Nsteps=Nsteps,Lmin=Lmin,mu=mu)

        xcoords[i,:] = np.array(xpos)
        ycoords[i,:] = np.array(ypos)
        Lss[i,:] = np.array(Ls)
        #print(f"Levies {i}")

    return (xcoords,ycoords,Lss)

def MAD(array):
    return np.median(np.abs(array - np.median(array)))


def LevyFlight(Nsets, Nruns, Nsteps, Lmin, mu, Tmin, nu, plot):
    Ntimes = Nsteps
    timePerStep = 1
    ts = np.linspace(0,Nsteps + Tmin * Nsteps * timePerStep, Ntimes)

    mads = np.zeros((Nsets,Ntimes))

    #bar = Bar(f"mu = {mu} Overall progress", max = Nsets)
    for i in tqdm(range(Nsets), desc="LevyFlight Progress", ):
        xss, yss, lss = getLevies(Nruns,Nsteps,Lmin,mu)


        #Convert step based index to time based index
        xts, yts = np.zeros((Nruns, Ntimes)), np.zeros((Nruns, Ntimes))

        for k in range(len(xss[:,0])):
            xs = xss[k,:]
            ys = yss[k,:]
            waittimes = np.array([randomT(Tmin=Tmin,nu=nu) for _ in range(len(xs))])
            current_time = 0

            xt, yt = np.array([None for _ in range(len(ts))]), np.array([None for _ in range(len(ts))])

            positionIndex = 0
            hasWaited = False
            x, y = 0, 0
            for j in range(len(ts)):
                while(xt[j] is None):
                    if current_time == ts[j]:
                        xt[j] = x
                        yt[j] = y
                        break

                    if current_time > ts[j] and hasWaited == True:
                        xt[j] = x
                        yt[j] = y
                        break

                    if current_time > ts[j] and hasWaited == False:
                        velocity = np.sqrt((xs[positionIndex] - xs[positionIndex-1])**2 + (ys[positionIndex] - ys[positionIndex-1])**2) / timePerStep
                        angle = np.arctan((ys[positionIndex] - ys[positionIndex-1]) / (xs[positionIndex] - xs[positionIndex-1]))
                        xt[j] = x - velocity * np.cos(angle) * (current_time-ts[j])
                        yt[j] = y - velocity * np.sin(angle) * (current_time-ts[j])                        
                        break

                    
                    if hasWaited == False:
                        current_time += waittimes[positionIndex]
                        hasWaited = True
                        continue

                    if hasWaited == True:
                        current_time += timePerStep
                        positionIndex += 1
                        ################################################################################## Kaj ko zmanjka pozicij?
                        x = xs[positionIndex]
                        y = ys[positionIndex]
                        hasWaited = False
                        continue
            
                        
            xts[k,:], yts[k,:] = xt, yt



        
        rss = np.sqrt(np.square(xts) + np.square(yts)) ##### Razdalja od izhodišča

        mad = []

        for j in range(len(rss[0,:])):
            mad.append(MAD(rss[:,j]))
            #print(f"Runs {j}")

        mads[i,:] = np.array(mad)

        #print(f"Sets {i}")
        #bar.next()

    #bar.finish()



    average_mad = np.average(mads, axis=0)
    stdev_mad =np.std(mads, axis = 0)
    ts = np.array(ts)

    variance2 = np.square(average_mad) / 0.45494
    variance2_error = 2 * average_mad * stdev_mad / 0.45494

    if plot:
        plt.plot(ts, variance2, color = "black", zorder = 2)
        plt.errorbar(ts,variance2, variance2_error, color = "gray", linestyle = "", zorder = 1)
        plt.xlabel("t")
        plt.ylabel("$\sigma^2$(t)")
        plt.show()

    lnVariance2 =  np.log(variance2[1:])
    lnVariance2_error = variance2_error[1:] / variance2[1:]
    lnts = np.log(ts[1:])

    

    x = sm.add_constant(lnts)
    model = sm.WLS(lnVariance2,x,weights=1/np.square(lnVariance2_error))
    results = model.fit()
    if plot:
        print(results.params)
        print(results.bse)

    if plot:
        plt.plot(lnts, lnVariance2, zorder=1)
        plt.errorbar(lnts,lnVariance2, lnVariance2_error, zorder=2)
        plt.plot(lnts, [results.params[0] + results.params[1]  * x for x in lnts], linestyle = "dashed", zorder=3, color = "black")
        plt.xlabel("ln(t)")
        plt.ylabel("ln($\sigma^2$(t))")
        plt.show()

    return (results.params[1], results.bse[1])


def LevyWalk(Nsets, Nruns, Nsteps, Lmin, mu, Tmin, nu, plot):
    velocity = 1
    Ntimes = Nsteps
    ts = np.linspace(0,Nsteps * Lmin / velocity + Nsteps * Tmin, Ntimes)

    mads = np.ones((Nsets,Ntimes))

    #bar = Bar("Overall progress", max = Nsets * Nruns) 
    for i in tqdm(range(Nsets), desc="LevyWalk Progress"):
        xss, yss, lss = getLevies(Nruns,Nsteps,Lmin,mu)

        #Convert step based index to time based index
        xts, yts = np.zeros((Nruns, Ntimes)), np.zeros((Nruns, Ntimes))

        for k in range(len(xss[:,0])):
            xs = xss[k,:]
            ys = yss[k,:]
            waittimes = np.array([randomT(Tmin=Tmin,nu=nu) for _ in range(len(xs))])
            current_time = 0

            xt, yt = np.array([None for _ in range(len(ts))]), np.array([None for _ in range(len(ts))])

            positionIndex = 0
            hasWaited = False
            x, y = 0, 0
            for j in range(len(ts)):
                while(xt[j] is None):
                    if current_time == ts[j]:
                        xt[j] = x
                        yt[j] = y
                        break

                    if current_time > ts[j] and hasWaited == True:
                        xt[j] = x
                        yt[j] = y
                        break

                    if current_time > ts[j] and hasWaited == False:
                        angle = np.arctan((ys[positionIndex] - ys[positionIndex-1]) / (xs[positionIndex] - xs[positionIndex-1]))
                        xt[j] = x - velocity * np.cos(angle) * (current_time-ts[j])
                        yt[j] = y - velocity * np.sin(angle) * (current_time-ts[j])                        
                        break

                    
                    if hasWaited == False:
                        current_time += waittimes[positionIndex]
                        hasWaited = True
                        continue

                    if hasWaited == True:
                        current_time += np.sqrt((xs[positionIndex+1] - xs[positionIndex])**2 + (ys[positionIndex+1] - ys[positionIndex])**2) / velocity
                        positionIndex += 1
                        ################################################################################## Kaj ko zmanjka pozicij?
                        x = xs[positionIndex]
                        y = ys[positionIndex]
                        hasWaited = False
                        continue
            
                        
            xts[k,:], yts[k,:] = xt, yt


        rss = np.sqrt(np.square(xts) + np.square(yts)) ##### Razdalja od izhodišča

        mad = []

        for j in range(len(rss[0,:])):
            mad.append(MAD(rss[:,j]))
            #print(f"Runs {j}")

        mads[i,:] = np.array(mad)

        #print(f"Sets {i}")

    #bar.finish()



    average_mad = np.average(mads, axis=0)
    stdev_mad =np.std(mads, axis = 0)
    ts = np.array(ts)

    variance2 = np.square(average_mad) / 0.45494
    variance2_error = 2 * average_mad * stdev_mad / 0.45494

    if plot:
        plt.plot(ts, variance2, color = "black", zorder = 2)
        plt.errorbar(ts,variance2, variance2_error, color = "gray", linestyle = "", zorder = 1)
        plt.xlabel("t")
        plt.ylabel("$\sigma^2$(t)")
        plt.show()
        

    lnVariance2 =  np.log(variance2[:])
    lnVariance2_error = variance2_error[:] / variance2[:]
    lnts = np.log(ts[:])

    #Remove woird values not compatible with logarithm
    badValues = (float("inf"), float("-inf"), float("NaN"))
    i = 0
    while i < len(lnVariance2):
        if lnVariance2[i] in badValues or lnVariance2_error[i] in badValues or lnts[i] in badValues:
            lnVariance2 = np.delete(lnVariance2, i)
            lnVariance2_error = np.delete(lnVariance2_error, i)
            lnts = np.delete(lnts, i)
            i -= 1
        i += 1



    x = sm.add_constant(lnts)
    model = sm.WLS(lnVariance2,x,weights=1/np.square(lnVariance2_error))
    results = model.fit()
    if plot:
        print(results.params)
        print(results.bse)

    if plot:
        plt.plot(lnts, lnVariance2, zorder=1)
        plt.errorbar(lnts,lnVariance2, lnVariance2_error, zorder=2)
        plt.plot(lnts, [results.params[0] + results.params[1]  * x for x in lnts], zorder=3, linestyle = "dashed", color = "black")
        plt.xlabel("ln(t)")
        plt.ylabel("ln($\sigma^2$(t))")
        plt.show()
    

    return (results.params[1], results.bse[1])


LevyFlight(Nsets=10, Nruns=10000, Nsteps=1000, Lmin=1, mu=2, Tmin=0.1, nu=1.5, plot = True)
LevyWalk(Nsets=10, Nruns=10000, Nsteps=1000, Lmin=1, mu=2, Tmin=0.1, nu=1.5, plot = True)


musflight = np.linspace(1.1, 10, 100)
flightnus = []
flighterror = []
for i in tqdm(range(len(musflight)), desc = "Overall Progress:"):

    val, error = LevyFlight(Nsets = 10, Nruns=1000, Nsteps=100, Lmin = 1, mu=musflight[i], plot=False)
    flightnus.append(val)
    flighterror.append(error)



mus = np.linspace(1.3, 10, 100)
walknus = []
walkerror = []
for i in tqdm(range(len(mus)), desc = "Overall Progress:"):
    
    val, error = LevyWalk(Nsets = 10, Nruns=100, Nsteps=1000, Lmin = 1, mu=mus[i], plot=False)
    walknus.append(val)
    walkerror.append(error)



plt.errorbar(musflight,flightnus,flighterror, label = "Simuliran $\\nu$($\mu$)")
mustheory = np.linspace(1.1, 10, 1000)
flighttheory = [2/(mu-1) if mu < 3 else 1 for mu in mustheory]
plt.plot(mustheory, flighttheory, linestyle = "dashed", label = "Teoretična vrednost")
plt.legend()
plt.ylabel("$\\nu$")
plt.xlabel("$\mu$")
plt.show()

flighttheory = np.array([2/(mu-1) if mu < 3 else 1 for mu in musflight])
plt.errorbar(musflight,np.abs(flightnus - flighttheory),flighterror, label = "Razlika med simulirano in teoretično vrednostjo")
plt.legend()
plt.ylabel("$\\nu$")
plt.xlabel("$\mu$")
plt.show()


plt.errorbar(mus,walknus,walkerror, label = "Simuliran $\\nu$($\mu$)")
mustheory = np.linspace(1.1, 10, 1000)
walktheory = [2 if mu < 2 else 4-mu if mu < 3 else 1 for mu in mustheory]
plt.plot(mustheory, walktheory, linestyle = "dashed", label = "Teoretična vrednost")
plt.legend()
plt.ylabel("$\\nu$")
plt.xlabel("$\mu$")
plt.show()

walktheory = np.array([2/(mu-1) if mu < 3 else 1 for mu in mus])
plt.errorbar(mus,np.abs(walknus - walktheory),walkerror, label = "Razlika med simulirano in teoretično vrednostjo")
plt.legend()
plt.ylabel("$\\nu$")
plt.xlabel("$\mu$")
plt.show()




fig, axs = plt.subplots(3,3)
i,j = 0,0
for mu in [1.1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5]:
    xdata, ydata, Ls = generateLevy(10000,1,mu)

    axs[i,j].plot(xdata,ydata, color = "gray")
    axs[i,j].scatter(xdata,ydata,s=10,color="black")
    axs[i,j].set_title(f"$\mu$ = {mu}")

    if j < 2:
        j+=1
    else:
        j=0
        i+=1

fig.subplots_adjust(left=0.035,
                    bottom=0.036, 
                    right=0.997, 
                    top=0.972, 
                    wspace=0.2, 
                    hspace=0.2)

plt.show()





    






