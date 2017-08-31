import math
import scipy.optimize

# Hydrogen frequency of the spectrometer
wa=850*10**6

# T_1/T_2
T1T2=11

gammaH=267.513*10**6
gammaN=-27.166*10**6
wx=wa*gammaN/gammaH
mu=4*3.14*10**(-7)
h_planck=1.055*10**(-34)
r=0.109*10**(-9)
rN=0.101*10**(-9)
d=(mu*gammaN*gammaH*h_planck)/(4*math.pi*rN**3)
c2=2*(math.pi*wx*160*10**(-6))**2/15

def J(w,t):
        return t/(1+(2*math.pi*w*t)**2) 

def f(t):
        return math.sqrt(((0.5*0.05*d**2*(4*J(0,t)+J(wa-wx,t)+3*J(wx,t)+6*J(wa,t)+6*J(wa+wx,t))+(2*math.pi*wx*160*10**(-6))**2*(3*J(wx,t)+4*J(0,t))/90)/(0.05*d**2*(J(wa-wx,t)+3*J(wx,t)+6*J(wa+wx,t))+(2*math.pi*wx*160*10**(-6))**2*J(wx,t)/15)-T1T2)**2)

sol = scipy.optimize.brute(f,[(10**(-10),1)])
print(sol)

        
