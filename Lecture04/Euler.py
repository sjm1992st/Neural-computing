# -*- coding: utf-8 -*-
from __future__ import division
from numpy import *
from pylab import *
import matplotlib.pylab as plt

class Neuron(object):
    def __init__(self,timeArray=None,inputCurrent=None):
        if timeArray is None:
            print "Using default time vector"

            timeArray = arange(0,850+mstep, mstep)
        if inputCurrent is None:
            print "Using default input current"
            inputCurrent = zeros(len(timeArray))
            for i,t in enumerate(timeArray):
                if t1 <= t <= t2:
                    inputCurrent[i] = Iinj

        self.add_v=[]
        self.timeArray = timeArray
        self.inputCurrent = inputCurrent

    def __call__(self):
        self.solveHHModel()
        self.plotVoltage()

        return self.voltageArray,self.dt,self.SpkNum

    def Euler_parameter(self,p_va,p_vb,kv):
         return p_va*(1-kv) - p_vb*kv
    def Euler_v(self,m,h,n,V,I):
        return (-gNa*m*m*m*h*(V-ENa)-gK*n*n*n*n*(V-EK)-gL*(V-EL)+I)/C

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

        m = 0.053
        h = 0.596
        n = 0.318

        for i in range(1,nt):
            I = inputCurrent[i-1]

            kv1=dt*self.Euler_v(m,h,n,V,I)
            km1=dt*self.Euler_parameter(alpham(V),betam(V),m)
            kh1=dt*self.Euler_parameter(alphah(V),betah(V),h)
            kn1=dt*self.Euler_parameter(alphan(V),betan(V),n)


            #print kv1,km1
            V = V+kv1
            m = m+km1
            h = h+kh1
            n = n+kn1
            #print V

            voltageArray[i] = V
            mSubunitArray[i] = m
            hSubunitArray[i] = h
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
            fig = plt.figure()
            ax = fig.add_subplot(111)
            ax.set_xlabel('Time (ms)')
            ax.set_ylabel('Voltage (mV)')
            ax.set_title('Euler_Firing rate: %d Hz' % SpkNum)
            plt.plot(t,V)
            show()

# Define subunit kinetics
def alpham(V):
    out = 0.32*(V+50)/(1-exp(-(V+50)*1.0/4))
    return out

def betam(V):
    out = -0.28*(V+23)/(1-exp((V+23)*1.0/5))
    return out

def alphah(V):
    out = 0.128*exp(-(V+46)*1.0/18)
    return out

def betah(V):
    out =4.0/(1+exp(-(V+23)*1.0/5))
    return out

def alphan(V):
    out = 0.032*(V+48)/(1-exp(-(V+48)*1.0/5))
    return out

def betan(V):
    out =0.5*exp(-(V+53)*1.0/40)
    return out

# End


#Define electrophysiological parameters
mstep=0.01
C = 0.2
gNa = 20
gK = 6
gL = 0.01
ENa = 50
EK = -90
EL = -65



##刺激时间
t1=100 #-ms
t2=520 #-ms


for Iinj in arange(1,20):      #Input current -uA
    hh = Neuron()
    print  Iinj
    mvoltageArray, mdt,mSpkNum=hh()
    print mSpkNum

