import sys
import numpy
import math
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

data1 = numpy.loadtxt(sys.argv[1], usecols=range(0,1)), numpy.loadtxt(sys.argv[1], usecols=range(1,2))
data2 = numpy.loadtxt(sys.argv[1], usecols=range(0,1)), numpy.loadtxt(sys.argv[1], usecols=range(2,3))
data3 = numpy.loadtxt(sys.argv[1], usecols=range(0,1)), numpy.loadtxt(sys.argv[1], usecols=range(3,4))

k1=[0]*3
k2=[0]*3
k3=[0]*3
kcov1=[0]*3
kcov2=[0]*3
kcov3=[0]*3
analL=2
for l in range(-1,2):
  j=0
  for i in range(0,len(data1[0])):
    if data1[0][i]<analL+l:
      j=j+1
      dataLOGx = [0]*j
      dataLOGy = [0]*j

  def func(x, k):
    return k*x
  #CFfig=plt.figure()

  for i in range(0,j):
    dataLOGx[i]=data1[0][i]
    dataLOGy[i]=data1[1][i]
  k1[l+1], kcov1[l+1] = curve_fit(func, dataLOGx, dataLOGy,p0=(1))
  #print k1,kcov1
  #plt.plot(dataLOGx,dataLOGy,dataLOGx,func(dataLOGx,k1))

  for i in range(0,j):
    dataLOGx[i]=data2[0][i]
    dataLOGy[i]=data2[1][i]
  k2[l+1], kcov2[l+1] = curve_fit(func, dataLOGx, dataLOGy,p0=(1))
  #print k2,kcov2
  #plt.plot(dataLOGx,dataLOGy,dataLOGx,func(dataLOGx,k2))

  for i in range(0,j):
    dataLOGx[i]=data3[0][i]
    dataLOGy[i]=data3[1][i]
  k3[l+1], kcov3[l+1] = curve_fit(func, dataLOGx, dataLOGy,p0=(1))
  #print k3,kcov3
  #plt.plot(dataLOGx,dataLOGy,dataLOGx,func(dataLOGx,k3))

#CFfig.savefig('../scratch/tst.pdf')

# RELAXATION TIMES ALONG PRINCIPAL AXES
#tau1=(-1)/k1
#tau2=(-1)/k2
#tau3=(-1)/k3
#print tau1,tau2,tau3

# DIFFUSION COEFFICIENTS
D1=(sum(k1)/len(k1))/2
D2=(sum(k2)/len(k2))/2
D3=(sum(k3)/len(k3))/2
# ERRORS FOR DIFFUSION CONSTANT
D1err=numpy.sqrt(((k1-D1*2) ** 2).mean())/2
D2err=numpy.sqrt(((k2-D2*2) ** 2).mean())/2
D3err=numpy.sqrt(((k3-D3*2) ** 2).mean())/2
RATIOerr=2*D1err/(D3+D2)+2*D1*D2err/(D3+D2)**2+2*D1*D3err/(D3+D2)**2

#print D1,D2,D3

# RELAXATION TIMES FOR N-H ROTATIONAL CORRELATION FUNCTION FROM OVERALL DIFFUSION 
tau1dot=1/(6*D3)
tau2dot=1/(D2+5*D3)
tau3dot=1/(4*D2+2*D3)

D=(D3+D2+D1)/3
L2=(D3*D2+D3*D1+D2*D3)/3
tau1=1/(4*D3+D2+D1)
tau2=1/(D3+4*D2+D1)
tau3=1/(D3+D2+4*D1)
if D**2-L2 < 0.0001:
  L2=D**2
tau4=1/(6*(D+math.sqrt(D**2-L2)))
tau5=1/(6*(D-math.sqrt(D**2-L2)))

#print tau1dot,tau2dot,tau3dot

print("D_xx        ", str(*(D3*100)), D3err*100)
print("D_yy        ", str(*(D2*100)), D2err*100)
print("D_zz        ", str(*(D1*100)), D1err*100)
print("D_||/D_+    ", str(*(2*D1/(D3+D2))), str(*(RATIOerr)))
print("D_av        ", str(*(D)*100),(D1err+D2err+D3err)*100/3)
print("tau1        ", str(*(tau1)))
print("tau2        ", str(*(tau2)))
print("tau3        ", str(*(tau3)))
print("tau4        ", str(*(tau4)))
print("tau5        ", str(*(tau5)))


#def fiveexpfunc(x,p1,p2,p3,p4,p5,tau1,tau2,tau3,tau4,tau5):
scalingF=1
def fiveexpfunc(x,p1,p2,p3,p4,p5):
  return p1**2*numpy.exp(-x/(scalingF*tau1)) + p2**2*numpy.exp(-x/(scalingF*tau2)) + p3**2*numpy.exp(-x/(scalingF*tau3)) + p4**2*numpy.exp(-x/(scalingF*tau4)) + p5**2*numpy.exp(-x/(scalingF*tau5))

for i in range(1,87):
  xdata = numpy.loadtxt(sys.argv[2]+'/overall/NHrotaCF_' + str(i) + '.xvg', usecols=range(0,1))
  ydata = numpy.loadtxt(sys.argv[2]+'/overall/NHrotaCF_' + str(i) + '.xvg', usecols=range(1,2))
  xdata = xdata*0.001

  j=0
  for k in range(0,len(xdata)):
    if xdata[k]<5:
      j=j+1
      dataLOGx = [0]*j
      dataLOGy = [0]*j

  for k in range(0,j):
    dataLOGx[k]=xdata[k]
    dataLOGy[k]=ydata[k]

  scalingF=1.0
  popt3, pcov3 = curve_fit(fiveexpfunc, dataLOGx, dataLOGy,p0=(0.45,0.45,0.45,0.45,0.45),maxfev=100000)
  print(popt3**2)

  #TSTfig=plt.figure()
  #axes = plt.gca()
  #axes.set_xlim([0,20])
  #plt.plot(xdata, ydata,xdata, twoexpfunc(xdata,*popt),xdata, threeexpfunc(xdata,*popt3))
  #plt.plot(xdata, ydata,xdata, fiveexpfunc(xdata,*popt3))
  #TSTfig.savefig('../scratch/tst2.pdf')

  xdata = numpy.loadtxt(sys.argv[2]+'/oriented/NHrotaCF_' + str(i) + '.xvg', usecols=range(0,1))
  ydata = numpy.loadtxt(sys.argv[2]+'/oriented/NHrotaCF_' +str(i)+'.xvg', usecols=range(1,2))
  xdata = xdata*0.001

  scalingF=float(sys.argv[3])
  ydata=ydata*fiveexpfunc(xdata,*popt3)

  #A=numpy.fft.fft([xdata,ydata])
  #print(A,2)
  
  scaledCFfile=open(str(sys.argv[2]+'/scaledrotation/NHrotaCF_' +str(i) + '.xvg'), 'w+')
  for j in range(0,len(xdata)):
    print(xdata[j]*1000,ydata[j], file=scaledCFfile)
                  
#TSTfig2=plt.figure()
#Axes = plt.gca()
#axes.set_xlim([0,20])
#plt.plot(xdata, ydata)
#TSTfig2.savefig('../scratch/tst3.pdf')

#def spectralDENS(x,p1,p2,p3,p4,p5,p6):
#  return 2*p1*p2/(2*2*math.pi*math.pi*p2*p2*x*x+1) + 2*p3*p4/(2*2*math.pi*math.pi*p4*p4*x*x+1) + 2*p5*p6/(2*2*math.pi*math.pi*p6*p6*x*x+1)

#def spectralDENS2(x,p1,p2,p3,p4,p5,p6):
#  return 2*p1*p2/(2*2*math.pi*math.pi*p2*p2*2*x*2*x+1) + 2*p3*p4/(2*2*math.pi*math.pi*p4*p4*2*x*2*x+1) + 2*p5*p6/(2*2*math.pi*math.pi*p6*p6*2*x*2*x+1)

 

#S=0.179
#tf=34*10**(-12)
#ts=3.8*10**(-9)

#def spectralDENSexp(x,S,tf,ts):
#  return 2*(1-S*S)*tf + 2*S*S*ts/(2*2*math.pi*math.pi*ts*ts*x*x+1)  

#def spectralDENSexp2(x,S,tf,ts):
#  return 2*(1-S*S)*tf + 2*S*S*ts/(2*2*math.pi*math.pi*ts*ts*2*x*2*x+1)  


#frequencyRANGE=numpy.arange(10**6, 10**8, 1000)

#ksi=167000
#def R1(x,p1,p2,p3,p4,p5,p6):
#  return 3*math.pi*math.pi*ksi*ksi*(2*spectralDENS(x,p1,p2,p3,p4,p5,p6)+8*spectralDENS2(x,p1,p2,p3,p4,p5,p6))/40

#def R1exp(x,S,tf,ts):
#  return 3*math.pi*math.pi*ksi*ksi*(2*spectralDENSexp(x,S,tf,ts)+8*spectralDENSexp2(x,S,tf,ts))/40

#def R2(x,p1,p2,p3,p4,p5,p6):
#  return 3*math.pi*math.pi*ksi*ksi*(3*spectralDENS(0,p1,p2,p3,p4,p5,p6)+5*spectralDENS(x,p1,p2,p3,p4,p5,p6)+2*spectralDENS2(x,p1,p2,p3,p4,p5,p6))/40

#def R2exp(x,S,tf,ts):
#  return 3*math.pi*math.pi*ksi*ksi*(3*spectralDENSexp(0,S,tf,ts)+5*spectralDENSexp(x,S,tf,ts)+2*spectralDENSexp2(x,S,tf,ts))/40


#Rfig=plt.figure()
#plt.plot(frequencyRANGE, R1exp(frequencyRANGE,S,tf,ts),frequencyRANGE, R1(frequencyRANGE,*popt3),Fdata,R1data)
#plt.plot(frequencyRANGE, R2exp(frequencyRANGE,S,tf,ts),frequencyRANGE, R2(frequencyRANGE,*popt3),Fdata,R2data)
#plt.xscale('log')
#Rfig.savefig('Rcomparison.pdf')

