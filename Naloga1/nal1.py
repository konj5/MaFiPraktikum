import numpy as np
import matplotlib.pyplot as plt
from mpmath import *

mp.dps = 500


#Taylor Series
def f(x:mpf, Nmax:int=1000) -> mpf:

    factor = 1
    sum = 0
    Nmin = 100

    for k in range(0,Nmax):
        #oklepaji k
        if k != 0:
            factor *= mpf(1/3+k-1)

        #3**k
        if k != 0:
            factor *= mpf(3)

        #x^3k
        if k != 0:
            factor *= x**mpf(3)

        #3k!
        if k != 0:
            factor /= mpf((3*k) * (3*k-1) * (3*k-2))

        new = sum + factor
        if new == sum and k > Nmin:
            break
        sum = new

    return sum

def g(x:mpf, Nmax:int=1000) -> mpf:

    factor = 1
    sum = 0
    Nmin = 100

    for k in range(0,Nmax):
        #oklepaji k
        if k != 0:
            factor *= mpf(2/3+k-1)
        #3**k
        if k != 0:
            factor *= mpf(3)

        #x^3k+1
        if k == 0:
            factor *= x
        else:
            factor *= x**mpf(3)

        #3k+1!
        if k != 0:
            factor /= mpf((3*k+1) * (3*k) * (3*k-1))

        new = sum + factor
        if new == sum and k > Nmin:
            break
        sum = new

    return sum

alfa = mpf(0.355028053887817239)
beta = mpf(0.258819403792806798)


#Asimptotic series

"""
def oldL(z:mpf, Nmax = 20) -> mpf:
    sum = 0
    factor = 1

    
    for s in range(0,Nmax):
        oldfactor = factor
        semifactor = 1
        #gamma(3s+1/2) / gamma(s+1/2)
        if s != 0:
            for k in range(s,3*s):
                semifactor *= mpf((k + 1/2))

        #s!
        if s != 0:
            factor /= mpf(s)

        #54^s
        if s != 0:
            factor /= mpf(54)

        #z^s
        if s != 0:
            factor /= mpf(z)

        #print(f"s={s}, {factor * semifactor}")

        #if oldfactor < factor * semifactor:break

        sum += factor * semifactor

    return sum

"""

def L(z:mpf, Nmax = 17) -> mpf:
    sum = 0
    factor = 1

    
    for s in range(0,Nmax):
        oldfactor = factor
        #gamma(3s+1/2) / gamma(s+1/2)
        if s != 0:
            factor *= mpf((3*s -3 + 1/2)) *  mpf((3*s -2 + 1/2)) *  mpf((3*s -1 + 1/2))

            factor /= mpf((s -1 + 1/2))

        #s!
        if s != 0:
            factor /= mpf(s)

        #54^s
        if s != 0:
            factor /= mpf(54)

        #z^s
        if s != 0:
            factor /= mpf(z)

        #print(f"s={s}, {factor}")

        #if oldfactor < factor:break

        sum += factor

    return sum

"""
def oldP(z:mpf, Nmax = 20) -> mpf:
    sum = 0
    factor = 1
    
    for s in range(0,Nmax):
        oldfactor = factor
        semifactor = 1
        #gamma(6s+1/2) / gamma(2s+1/2)
        if s != 0:
            for k in range(2*s,6*s):
                semifactor *= mpf((k + 1/2))

        #(2s)!
        if s != 0:
            factor /= mpf(2*s) * mpf(2*s-1)

        #54^2s
        if s != 0:
            factor /= mpf(54)**mpf(2)

        #z^2s
        if s != 0:
            factor /= mpf(z)**mpf(2)

        #(-1)^s
        if s != 0:
            factor *= mpf(-1)

        #print(f"s={s}, {factor * semifactor}")

        #if oldfactor < factor * semifactor:break

        sum += factor * semifactor

    return sum
"""

def P(z:mpf, Nmax = 23) -> mpf:
    sum = 0
    factor = 1
    
    for s in range(0,Nmax):
        oldfactor = factor
        #gamma(6s+1/2) / gamma(2s+1/2)
        if s != 0:
            for k in range(6*(s-1), 6*s):
                factor *= mpf((k + 1/2))
            factor /= (2*s-2+1/2) * (2*s-1+1/2)

        #(2s)!
        if s != 0:
            factor /= mpf(2*s) * mpf(2*s-1)

        #54^2s
        if s != 0:
            factor /= mpf(54)**mpf(2)

        #z^2s
        if s != 0:
            factor /= mpf(z)**mpf(2)

        #(-1)^s
        if s != 0:
            factor *= mpf(-1)

        #print(f"s={s}, {factor}")

        #if oldfactor < factor * semifactor:break

        sum += factor

    return sum

"""
def oldQ(z:mpf, Nmax = 3) -> mpf:
    sum = 0
    factor = 1
    
    for s in range(0,Nmax):
        oldfactor = factor
        semifactor = 1
        #gamma(3(2s+1)+1/2) / gamma((2s+1)+1/2)
        for k in range(2*s+1,6*s+3):
            semifactor *= mpf((k + 1/2))

        #(2s+1)!
        if s != 0:
            factor /= mpf(2*s+1) * mpf(2*s)

        #54^2s+1
        if s == 0:
            factor /= mpf(54)
        else:
            factor /= mpf(54)**mpf(2)

        #z^2s+1
        if s == 0:
            factor /= mpf(z)
        else:
            factor /= mpf(z)**mpf(2)

        #(-1)^s
        if s != 0:
            factor *= mpf(-1)

        #print(f"s={s}, {factor * semifactor}")

        #if oldfactor < factor * semifactor:break

        sum += factor * semifactor

    return sum
"""

def Q(z:mpf, Nmax = 20) -> mpf:
    sum = 0
    factor = 1
    
    for s in range(0,Nmax):
        oldfactor = factor
        #gamma(3(2s+1)+1/2) / gamma((2s+1)+1/2)
        if s != 0:
            for k in range(-3, 3):
                factor *= mpf((6*s + k + 1/2))
            factor /= (2*s-1+1/2) * (2*s + 1/2)
        if s == 0:
            factor *= (1+1/2)*(2+1/2)

        #(2s+1)!
        if s != 0:
            factor /= mpf(2*s+1) * mpf(2*s)

        #54^2s+1
        if s == 0:
            factor /= mpf(54)
        else:
            factor /= mpf(54)**mpf(2)

        #z^2s+1
        if s == 0:
            factor /= mpf(z)
        else:
            factor /= mpf(z)**mpf(2)

        #(-1)^s
        if s != 0:
            factor *= mpf(-1)

        #print(f"s={s}, {factor}")

        #if oldfactor < factor * semifactor:break

        sum += factor

    return sum
 


def Ai(x:mpf, type:str = "compound") -> mpf:
    if type == "taylor":
        return alfa * f(x) - beta * g(x)
    
    if type == "asimp":
        ksi = mpf(2/3)*abs(x)**mpf(3/2)

        if x > 0: 
            return mpf(np.e)**(-ksi) / (mpf(2) * mpf(np.pi)**mpf(1/2) * x**mpf(1/4)) * L(-ksi)
        
        if x < 0:
            return 1/(mpf(np.pi)**mpf(1/2) * (-x)**mpf(1/4)) * (sin(ksi - mpf(np.pi/4)) * Q(ksi) + cos(ksi - mpf(np.pi/4)) * P(ksi))
    
    if type == "compound":
        if x < -7.85 or x > 5.54:
            return Ai(x,"asimp")
        else:
            return Ai(x,"taylor")
            
    
def Bi(x:mpf, type:str = "compound") -> mpf:
    if type == "taylor":
        return np.sqrt(3) * (alfa * f(x) + beta * g(x))
    
    if type == "asimp":
        ksi = mpf(2/3)*abs(x)**mpf(3/2)

        if x > 0: 
            return mpf(np.e)**(ksi) / (mpf(np.pi)**mpf(1/2) * x**mpf(1/4)) * L(ksi)
        
        if x < 0:
            return 1/(mpf(np.pi)**mpf(1/2) * (-x)**mpf(1/4)) * (-sin(ksi - mpf(np.pi/4)) * P(ksi) + cos(ksi - mpf(np.pi/4)) * Q(ksi))

    if type == "compound":
        if x < -7.3 or x > 100:
            return Bi(x,"asimp")
        else:
            return Bi(x,"taylor")


xs = np.linspace(-80,80,1000)


AiTaylorAbsError = np.array([abs(airyai(mpf(x)) - Ai(mpf(x), "taylor")) for x in xs])
BiTaylorAbsError = np.array([abs(airybi(mpf(x)) - Bi(mpf(x), "taylor")) for x in xs])

AiAsimpAbsError = np.array([abs(airyai(mpf(x)) - Ai(mpf(x), "asimp")) for x in xs])
BiAsimpAbsError = np.array([abs(airybi(mpf(x)) - Bi(mpf(x), "asimp")) for x in xs])


plt.plot(xs, AiTaylorAbsError, label = "Ai - $Ai_{Taylor}$")
plt.plot(xs, AiAsimpAbsError,  label = "Ai - $Ai_{Asimp}$")
plt.axhline(10**-10, linestyle = "dashed", color = "black")
plt.legend()
plt.yscale("log")
plt.xlabel("x")
plt.grid()
plt.show()

plt.plot(xs, BiTaylorAbsError, label = "Bi - $Bi_{Taylor}$")
plt.plot(xs, BiAsimpAbsError,  label = "Bi - $Bi_{Asimp}$")
plt.axhline(10**-10, linestyle = "dashed", color = "black")
plt.legend()
plt.yscale("log")
plt.xlabel("x")
plt.grid()
plt.show()



AiAbsError = np.array([abs(airyai(mpf(x)) - Ai(mpf(x), "compound")) for x in xs])
BiAbsError = np.array([abs(airybi(mpf(x)) - Bi(mpf(x), "compound")) for x in xs])

AiRelError = np.array([abs(1 - Ai(mpf(x), "compound")/airyai(mpf(x))) for x in xs])
BiRelError = np.array([abs(1- Bi(mpf(x), "compound")/airybi(mpf(x))) for x in xs])


plt.plot(xs, AiAbsError, label = "Absolutna napaka")
plt.plot(xs, AiRelError,  label = "Relativna napaka")
plt.axhline(10**-10, linestyle = "dashed", color = "black")
plt.legend()
plt.xlabel("x")
plt.yscale("log")
plt.grid()
plt.show()

plt.plot(xs, BiAbsError, label = "Absolutna napaka")
plt.plot(xs, BiRelError,  label = "Relativna napaka")
plt.axhline(10**-10, linestyle = "dashed", color = "black")
plt.legend()
plt.xlabel("x")
plt.yscale("log")
plt.grid()
plt.show()


def bisection(f:callable, xmin:float, xmax:float, tolerance:float = mpf(10)**-10):
    assert sign(f(xmin)) * sign(f(xmax)) < 0

    xmid = (xmax + xmin)/2
    if abs(xmax - xmin) < tolerance: return xmid

    midval = f(xmid)
    

    if sign(f(xmin)) * sign(midval) < 0:
        return bisection(f,xmin,xmid)
    
    if sign(f(xmax)) * sign(midval) < 0:
        return bisection(f,xmid,xmax)
    
def func(z:float):
    return z**(2/3) * (1 + 5/48 / z**2 - 5/36 / z**4 + 77125/82944 / z**6 - 108056875/6967296 / z**8)


for stuff in [(airyai, 1), (airybi,3)]:
    zeros = []
    step = mpf(0.1)
    xmin = mpf(0)
    xmax = -step
    f = stuff[0]

    zerosapprox = []
    for s in range(1,101):
        zerosapprox.append(-func(3 * np.pi / 8 * (4*s - stuff[1])))

    while(len(zeros) < 100):
        print(xmin)
        if sign(f(xmin)) * sign(f(xmax)) < 0:
            zeros.append(bisection(f,xmin,xmax))
        
        xmin -= step
        xmax -= step

    #print(zeros)


    xs = np.linspace(0, float(zeros[-1]) - 1, 1000)
    fs = [f(x) for x in xs]
    plt.plot(xs,fs)
    plt.scatter(zeros, [0 for _ in range(len(zeros))], color = "red")
    plt.grid()
    plt.xlabel("x")
    plt.show()

    ns = [i for i in range(1,101)]
    error = [abs(zeros[i]-zerosapprox[i]) for i in range(len(ns))]
    plt.scatter(ns, error)
    plt.grid()
    plt.xlabel("Zaporedna številka ničle")
    plt.yscale("log")
    plt.show()
