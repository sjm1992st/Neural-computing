# -*- coding: utf-8 -*-
from __future__ import division
from numpy import *
from pylab import *
import matplotlib.pylab as plt
import numpy.linalg as nplg
class Neuron(object):
    def __init__(self,timeArray=None,inputCurrent=None):
        if timeArray is None:
            print "Using default time vector"

            timeArray = arange(0,300+mstep, mstep)
        if inputCurrent is None:
            print "Using default input current"
            inputCurrent = zeros(len(timeArray))
            for i,t in enumerate(timeArray):
                if t1 <= t <= t2:
                    inputCurrent[i] = t*1.0/30

        self.add_v=[]
        self.timeArray = timeArray
        self.inputCurrent = inputCurrent

    def __call__(self):
        self.solveHHModel()
        self.plotVoltage()

        return self.voltageArray,self.dt,self.SpkNum

    def Euler_v(self,m,n,V,I):
        return (-gNa*m*(V-ENa)-gK*n*(V-EK)-gL*(V-EL)+I)/C

    def solveHHModel(self):
        t = self.timeArray
        inputCurrent = self.inputCurrent

        dt = t[1] - t[0]
        nt = len(t)
        self.dt=dt
        #Initialize output arrays
        voltageArray = zeros((nt))
        mSubunitArray = zeros((nt))
        hSubunitArray = zeros((nt))
        nSubunitArray = zeros((nt))


        t0 = t[0]
        V = -65

        n = 0
        for i in range(1,nt):
            I = inputCurrent[i-1]

            kv1=dt*self.Euler_v(m_inf(V),n,V,I)
            kn1=dt*(n_inf(V)-n)
            kv2=dt*self.Euler_v(m_inf(V+kv1),n+kn1,V+kv1,I)
            kn2=dt*(n_inf(V+kv1)-n)


            #print kv1,km1
            V = V+(kv1+kv2)/2
            n = n+(kn1+kn2)/2
            #print V

            voltageArray[i] = V
            nSubunitArray[i] = n

        self.voltageArray = voltageArray
        self.mSubunitArray = mSubunitArray
        self.hSubunitArray = hSubunitArray
        self.nSubunitArray = nSubunitArray


    def plotVoltage(self,ax = None,lineStyle = 'k'):
        t = self.timeArray
        V = self.voltageArray
        sti = self.inputCurrent

        count = 0
        for i in range(int(t1/mstep),int(t2/mstep)):
            if (V[i]>V[i-1]) and (V[i]>V[i+1]) and V[i]>0:
                count += 1
        ##print count
        SpkNum = count*1.0*1000/(t2-t1)
        self.SpkNum = SpkNum


        if ax is None:
            plt.subplot(211)
            plt.xlabel('Time (ms)')
            plt.ylabel('Voltage (mV)')
            plt.title('Heun_Firing rate: %d Hz' % SpkNum)
            plt.plot(t,V)
            plt.subplot(212)
            plt.xlabel('Time (ms)')
            plt.ylabel('Injected current(uA)')
            plt.plot(t,sti)
            show()

# Define subunit kinetics
def m_inf(V):
    out = 1.0/(1+exp((-20-V)*1.0/15))
    return out

def n_inf(V):
    out = 1.0/(1+exp((-30-V)*1.0/5))
    return out


# End


#Define electrophysiological parameters
mstep=0.01
C = 1
gNa = 20
gK = 10
gL = 8
ENa = 60
EK = -90
EL = -80


##刺激时间
t1=0
t2=300
def trv(V,n):
    Ia=1
    dvdt=(-gNa*m_inf(V)*(V-ENa)-gK*n*(V-EK)-gL*(V-EL)+Ia)/C
    print dvdt
    dndt=n_inf(V)-n
    a=-gL-gNa*m_inf(V)-gNa*(V-ENa)*(exp((-20-V)*1.0/15)*1.0/15)*1.0*(m_inf(V)**2)-gK*n-gK*(V-EK)*(n_inf(V)-n)*1.0/dvdt
    #print a,dvdt,dndt
    b=a*dvdt*1.0/dndt
    c=0.2*exp((-30-V)*1.0/5)*(n_inf(V)**2)-dndt*1.0/(dvdt)
    d=0.2*exp((-30-V)*1.0/5)*(n_inf(V)**2)*dndt*1.0/(dvdt)-1
    print a,b,c,d
    trz=np.array([[a,b],[c,d]])
    #print trz.shape
    print nplg.eig(trz)
    print (a+d)**2-4*(a*d-b*c)

hh = Neuron()
#print np.dot(trz,np.array([[0.77519965],[0.63171632]])),3.25252725e+08*np.array([[0.77519965],[0.63171632]])
mvoltageArray, mdt,mSpkNum=hh()
print mSpkNum

V = -65.3046
n=0.00135667
trv(V,n)

V = -55.9619
n=0.0056563
trv(V,n)

V = -34.68
n=0.2813
trv(V,n)
#plt.annotate('',xy=(2,2),xytext=(8,8),arrowprops=dict(arrowstyle="<-",connectionstyle="arc3"))
##########
#########
