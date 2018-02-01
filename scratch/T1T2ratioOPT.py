import math
import scipy.optimize
import sys

# Hydrogen frequency of the spectrometer
wa=850*10**6

# T_1/T_2
#T1T2=float(sys.argv[1])

gammaH=267.513*10**6
gammaN=-27.166*10**6
wx=wa*gammaN/gammaH
mu=4*3.14*10**(-7)
h_planck=1.055*10**(-34)
r=0.109*10**(-9)
rN=0.101*10**(-9)
d=(mu*gammaN*gammaH*h_planck)/(4*math.pi*rN**3)
c2=2*(math.pi*wx*160*10**(-6))**2/15

solSUM=0
solSTDsum=0
sum1=0

def J(w,t):
        return t/(1+(2*math.pi*w*t)**2) 

crs=open(sys.argv[1],"r")
for columns in ( raw.strip().split() for raw in crs ):  
    T1T2=float(columns[1])
    def f(t):
        return math.sqrt(((0.5*0.05*d**2*(4*J(0,t)+J(wa-wx,t)+3*J(wx,t)+6*J(wa,t)+6*J(wa+wx,t))+(2*math.pi*wx*160*10**(-6))**2*(3*J(wx,t)+4*J(0,t))/90)/(0.05*d**2*(J(wa-wx,t)+3*J(wx,t)+6*J(wa+wx,t))+(2*math.pi*wx*160*10**(-6))**2*J(wx,t)/15)-T1T2)**2)
    sum1=sum1+1
    sol = scipy.optimize.brute(f,[(10**(-10),1)])
    #print(float(columns[0]),float(sol))
    solSUM=solSUM+sol
crs.close()
solAV=solSUM/sum1
crs=open(sys.argv[1],"r")
for columns in ( raw.strip().split() for raw in crs ):  
    T1T2=float(columns[1])
    def f(t):
        return math.sqrt(((0.5*0.05*d**2*(4*J(0,t)+J(wa-wx,t)+3*J(wx,t)+6*J(wa,t)+6*J(wa+wx,t))+(2*math.pi*wx*160*10**(-6))**2*(3*J(wx,t)+4*J(0,t))/90)/(0.05*d**2*(J(wa-wx,t)+3*J(wx,t)+6*J(wa+wx,t))+(2*math.pi*wx*160*10**(-6))**2*J(wx,t)/15)-T1T2)**2)
    sol = scipy.optimize.brute(f,[(10**(-10),1)])
    #print((sol-solAV)**2)
    solSTDsum=solSTDsum+(sol-solAV)**2
crs.close()

solSTD=math.sqrt(solSTDsum/(sum1-1))
print(solAV,solSTD/math.sqrt(sum1))
        
