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

    def RungeKutta_parameter(self,p_va,p_vb,kv):
         return p_va*(1-kv) - p_vb*kv
    def RungeKutta_v(self,m,h,n,V,I):
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

            kv1=dt*self.RungeKutta_v(m,h,n,V,I)
            km1=dt*self.RungeKutta_parameter(alpham(V),betam(V),m)
            kh1=dt*self.RungeKutta_parameter(alphah(V),betah(V),h)
            kn1=dt*self.RungeKutta_parameter(alphan(V),betan(V),n)

            kv2=dt*self.RungeKutta_v(m+0.5*km1,h+0.5*kh1,n+0.5*kn1,V+0.5*kv1,I)
            km2=dt*self.RungeKutta_parameter(alpham(V+0.5*kv1),betam(V+0.5*kv1),m+0.5*km1)
            kh2=dt*self.RungeKutta_parameter(alphah(V+0.5*kv1),betah(V+0.5*kv1),h+0.5*kh1)
            kn2=dt*self.RungeKutta_parameter(alphan(V+0.5*kv1),betan(V+0.5*kv1),n+0.5*kn1)

            kv3=dt*self.RungeKutta_v(m+0.5*km2,h+0.5*kh2,n+0.5*kn2,V+0.5*kv2,I)
            km3=dt*self.RungeKutta_parameter(alpham(V+0.5*kv2),betam(V+0.5*kv2),m+0.5*km2)
            kh3=dt*self.RungeKutta_parameter(alphah(V+0.5*kv2),betah(V+0.5*kv2),h+0.5*kh2)
            kn3=dt*self.RungeKutta_parameter(alphan(V+0.5*kv2),betan(V+0.5*kv2),n+0.5*kn2)

            kv4=dt*self.RungeKutta_v(m+km3,h+kh3,n+kn3,V+kv3,I)
            km4=dt*self.RungeKutta_parameter(alpham(V+kv3),betam(V+kv3),m+km3)
            kh4=dt*self.RungeKutta_parameter(alphah(V+kv3),betah(V+kv3),h+kh3)
            kn4=dt*self.RungeKutta_parameter(alphan(V+kv3),betan(V+kv3),n+kn3)

            V = V+(kv1+2*kv2+2*kv3+kv4)/6
            m = m+(km1+2*km2+2*km3+km4)/6
            h = h+(kh1+2*kh2+2*kh3+kh4)/6
            n = n+(kn1+2*kn2+2*kn3+kn4)/6
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
            ax.set_title('RK_Firing rate: %d Hz' % SpkNum)
            plt.plot(t,V)
            show()

# Define subunit kinetics
def alpham(V):
    out = 0.1*(V+40)/(1-exp(-(V+40)*1.0/10))
    return out

def betam(V):
    out = 4*exp(-(V+65)*1.0/18)
    return out

def alphah(V):
    out = 0.07*exp(-(V+65)*1.0/20)
    return out

def betah(V):
    out = 1.0/(1+exp(-(V+35)*1.0/10))
    return out

def alphan(V):
    out = 0.01*(V+55)/(1-exp(-(V+55)*1.0/10))
    return out

def betan(V):
    out = 0.125*exp(-(V+65)*1.0/80)
    return out

# End


#Define electrophysiological parameters
mstep=0.01
C = 1
gNa = 120
gK = 36
gL = 0.3
ENa = 50
EK = -77
EL = -54.4



##刺激时间
t1=100
t2=520


for Iinj in arange(7,62.1,5.51):      #Input current
    hh = Neuron()
    print  Iinj
    mvoltageArray, mdt,mSpkNum=hh()
    print mSpkNum





