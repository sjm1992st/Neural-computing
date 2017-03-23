# -*- coding: utf-8 -*-
from __future__ import division
from numpy import *
from pylab import *
import matplotlib.pylab as plt

class Neuron(object):
    def __init__(self,timeArray=None,inputCurrent=None):
        if timeArray is None:
            print "Using default time vector"

            timeArray = arange(0,500+mstep, mstep)
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
    def Euler_v(self,m,h,n,p,q,V,I):
        return (-gNa*m*m*m*h*(V-ENa)-gK*n*n*n*n*(V-EK)-gCa*p*p*p*q*(V-ECa)-gL*(V-EL)+I)/C

    def solveHHModel(self):
        t = self.timeArray
        inputCurrent = self.inputCurrent

        dt = t[1] - t[0]
        nt = len(t)
        self.dt=dt
        #Initialize output arrays
        voltageArray = zeros((nt))
        amSubunitArray = zeros((nt))
        ahSubunitArray = zeros((nt))
        anSubunitArray = zeros((nt))
        apSubunitArray = zeros((nt))
        aqSubunitArray = zeros((nt))
        bmSubunitArray = zeros((nt))
        bhSubunitArray = zeros((nt))
        bnSubunitArray = zeros((nt))
        bpSubunitArray = zeros((nt))
        bqSubunitArray = zeros((nt))

        t0 = t[0]
        V = -65

        m = 0
        h = 0
        n = 0
        p=0
        q=0

        for i in range(1,nt):
            I = inputCurrent[i-1]

            voltageArray[i] = V
            amSubunitArray[i] = alpham(V)
            ahSubunitArray[i] = alphah(V)
            anSubunitArray[i] = alphan(V)
            apSubunitArray[i] = alphap(V)
            aqSubunitArray[i] = alphaq(V)
            bmSubunitArray[i] = betam(V)
            bhSubunitArray[i] = betah(V)
            bnSubunitArray[i] = betan(V)
            bpSubunitArray[i] = betap(V)
            bqSubunitArray[i] = betaq(V)

            kv1=dt*self.Euler_v(m,h,n,p,q,V,I)
            km1=dt*self.Euler_parameter(alpham(V),betam(V),m)
            kh1=dt*self.Euler_parameter(alphah(V),betah(V),h)
            kn1=dt*self.Euler_parameter(alphan(V),betan(V),n)
            kp1=dt*self.Euler_parameter(alphap(V),betap(V),p)
            kq1=dt*self.Euler_parameter(alphaq(V),betaq(V),q)

            #print kv1,km1
            kv2=dt*self.Euler_v(m+km1,h+kh1,n+kn1,p+kp1,q+kq1,V+kv1,I)
            km2=dt*self.Euler_parameter(alpham(V+kv1),betam(V+kv1),m+km1)
            kh2=dt*self.Euler_parameter(alphah(V+kv1),betah(V+kv1),h+kh1)
            kn2=dt*self.Euler_parameter(alphan(V+kv1),betan(V+kv1),n+kn1)
            kp2=dt*self.Euler_parameter(alphah(V+kv1),betah(V+kv1),p+kp1)
            kq2=dt*self.Euler_parameter(alphan(V+kv1),betan(V+kv1),q+kq1)

            #print kv1,km1
            V = V+(kv1+kv2)/2
            m = m+(km1+km2)/2
            h = h+(kh1+kh2)/2
            n = n+(kn1+kn2)/2
            p = p+(kp1+kp2)/2
            q = q+(kq1+kq2)/2

            """
            V = V+kv1
            m = m+km1
            h = h+kh1
            n = n+kn1
            p = p+kp1
            q = q+kq1
            """
            #print V
           # voltageArray[i] = V


        self.voltageArray = voltageArray
        self.amSubunitArray = amSubunitArray
        self.ahSubunitArray = ahSubunitArray
        self.anSubunitArray = anSubunitArray
        self.apSubunitArray = apSubunitArray
        self.aqSubunitArray = aqSubunitArray
        self.bmSubunitArray = bmSubunitArray
        self.bhSubunitArray = bhSubunitArray
        self.bnSubunitArray = bnSubunitArray
        self.bpSubunitArray = bpSubunitArray
        self.bqSubunitArray = bqSubunitArray

    def plotVoltage(self,ax = None,lineStyle = 'k'):
        t = self.timeArray
        V = self.voltageArray
        sti = self.inputCurrent
        am = self.amSubunitArray
        an = self.anSubunitArray
        ah = self.ahSubunitArray
        ap = self.apSubunitArray
        aq = self.aqSubunitArray
        bm = self.bmSubunitArray
        bn = self.bnSubunitArray
        bh = self.bhSubunitArray
        bp = self.bpSubunitArray
        bq = self.bqSubunitArray

        count = 0
        for i in range(int(t1/mstep),int(t2/mstep)):
            if (V[i]>V[i-1]) and (V[i]>V[i+1]) and V[i]>0:
                count += 1
        ##print count
        SpkNum = count*1.0*1000/(t2-t1)
        self.SpkNum = SpkNum




        kflag=zeros((2))
        flagcount=0
        for j in  range(len(V)-2):
            if V[j+1]<V[j] and V[j+1]<V[j+2]:
                kflag[flagcount]=j+1
                flagcount=flagcount+1
                #print V[kflag-1],V[kflag],V[kflag+1]

            if flagcount==2:
                break
        if ax is None:
            fig = plt.figure()
            ax = fig.add_subplot(111)
            ax.set_xlabel('Time (ms)')
            ax.set_ylabel('Voltage (mV)')
            ax.set_title('Euler_Firing rate: %d Hz' % SpkNum)
            plt.plot(t,V)

            fig = plt.figure()
            ax = fig.add_subplot(5,2,1)
            ax.set_xlabel('Voltage (mV)')
            ax.set_ylabel('Rate')
            ax.set_title('alpha-beta')
            plt.plot(V[kflag[0]:kflag[1]],am[kflag[0]:kflag[1]],'b',label='alpham')
            plt.plot(V[kflag[0]:kflag[1]],bm[kflag[0]:kflag[1]],'g',label='betam')
            plt.annotate('am', xy = (V[40],am[40]), xytext = (V[40],am[40]))
            plt.annotate('bm', xy = (V[40],bm[40]), xytext = (V[40],bm[40]))
            #print V[kflag[0]:kflag[1]]
            ax = fig.add_subplot(5,2,2)
            ax.set_xlabel('Voltage (mV)')
            ax.set_ylabel('Time constant')
            ax.set_title('tao-inf')
            tao=1.0/(am[kflag[0]:kflag[1]]+bm[kflag[0]:kflag[1]])
            inf=am[kflag[0]:kflag[1]]*1.0/(am[kflag[0]:kflag[1]]+bm[kflag[0]:kflag[1]])
            ax.plot(V[kflag[0]:kflag[1]],tao,'b--',label='tao')
            ax2=ax.twinx()
            ax2.plot(V[kflag[0]:kflag[1]],inf,'g',label='inf')
            ax.annotate('$\\tau$', xy = (V[40],tao[40]), xytext = (V[40],tao[40]))
            ax2.annotate('m$\infty$', xy = (V[40],inf[40]), xytext = (V[40],inf[40]))

            ax = fig.add_subplot(5,2,3)
            plt.plot(V[kflag[0]:kflag[1]],an[kflag[0]:kflag[1]],'b',label='alphan')
            plt.plot(V[kflag[0]:kflag[1]],bn[kflag[0]:kflag[1]],'g',label='betan')
            plt.annotate('an', xy = (V[40],an[40]), xytext = (V[40],an[40]))
            plt.annotate('bn', xy = (V[40],bn[40]), xytext = (V[40],bn[40]))
            #print V[kflag[0]:kflag[1]]
            ax = fig.add_subplot(5,2,4)
            tao=1.0/(an[kflag[0]:kflag[1]]+bn[kflag[0]:kflag[1]])
            inf=an[kflag[0]:kflag[1]]*1.0/(an[kflag[0]:kflag[1]]+bn[kflag[0]:kflag[1]])
            ax.plot(V[kflag[0]:kflag[1]],tao,'b--',label='tao')
            ax2=ax.twinx()
            ax2.plot(V[kflag[0]:kflag[1]],inf,'g',label='inf')
            ax.annotate('$\\tau$', xy = (V[40],tao[40]), xytext = (V[40],tao[40]))
            ax2.annotate('n$\infty$', xy = (V[40],inf[40]), xytext = (V[40],inf[40]))

            ax = fig.add_subplot(5,2,5)
            plt.plot(V[kflag[0]:kflag[1]],ah[kflag[0]:kflag[1]],'b',label='alphah')
            plt.plot(V[kflag[0]:kflag[1]],bh[kflag[0]:kflag[1]],'g',label='betah')
            plt.annotate('ah', xy = (V[40],ah[40]), xytext = (V[40],ah[40]))
            plt.annotate('bh', xy = (V[40],bh[40]), xytext = (V[40],bh[40]))
            #print V[kflag[0]:kflag[1]]
            ax = fig.add_subplot(5,2,6)
            tao=1.0/(ah[kflag[0]:kflag[1]]+bh[kflag[0]:kflag[1]])
            inf=ah[kflag[0]:kflag[1]]*1.0/(ah[kflag[0]:kflag[1]]+bh[kflag[0]:kflag[1]])
            ax.plot(V[kflag[0]:kflag[1]],tao,'b--',label='tao')
            ax2=ax.twinx()
            ax2.plot(V[kflag[0]:kflag[1]],inf,'g',label='inf')
            ax.annotate('$\\tau$', xy = (V[40],tao[40]), xytext = (V[40],tao[40]))
            ax2.annotate('h$\infty$', xy = (V[40],inf[40]), xytext = (V[40],inf[40]))

            ax = fig.add_subplot(5,2,7)
            plt.plot(V[kflag[0]:kflag[1]],ap[kflag[0]:kflag[1]],'b',label='alphap')
            plt.plot(V[kflag[0]:kflag[1]],bp[kflag[0]:kflag[1]],'g',label='betap')
            plt.annotate('ap', xy = (V[40],ap[40]), xytext = (V[40],ap[40]))
            plt.annotate('bp', xy = (V[40],bp[40]), xytext = (V[40],bp[40]))
            #print V[kflag[0]:kflag[1]]
            ax = fig.add_subplot(5,2,8)
            tao=1.0/(ap[kflag[0]:kflag[1]]+bp[kflag[0]:kflag[1]])
            inf=ap[kflag[0]:kflag[1]]*1.0/(ap[kflag[0]:kflag[1]]+bp[kflag[0]:kflag[1]])
            ax.plot(V[kflag[0]:kflag[1]],tao,'b--',label='tao')
            ax2=ax.twinx()
            ax2.plot(V[kflag[0]:kflag[1]],inf,'g',label='inf')
            ax.annotate('$\\tau$', xy = (V[40],tao[40]), xytext = (V[40],tao[40]))
            ax2.annotate('p$\infty$', xy = (V[40],inf[40]), xytext = (V[40],inf[40]))

            ax = fig.add_subplot(5,2,9)
            plt.plot(V[kflag[0]:kflag[1]],aq[kflag[0]:kflag[1]],'b',label='alphaq')
            plt.plot(V[kflag[0]:kflag[1]],bq[kflag[0]:kflag[1]],'g',label='betaq')
            ax.set_xlabel('Voltage (mV)')
            plt.annotate('aq', xy = (V[40],aq[40]), xytext = (V[40],aq[40]))
            plt.annotate('bq', xy = (V[40],bq[40]), xytext = (V[40],bq[40]))
            #print V[kflag[0]:kflag[1]]
            ax = fig.add_subplot(5,2,10)
            tao=1.0/(aq[kflag[0]:kflag[1]]+bq[kflag[0]:kflag[1]])
            inf=aq[kflag[0]:kflag[1]]*1.0/(aq[kflag[0]:kflag[1]]+bq[kflag[0]:kflag[1]])
            ax.plot(V[kflag[0]:kflag[1]],tao,'b--',label='tao')
            ax2=ax.twinx()
            ax2.plot(V[kflag[0]:kflag[1]],inf,'g',label='inf')
            ax.annotate('$\\tau$', xy = (V[40],tao[40]), xytext = (V[40],tao[40]))
            ax2.annotate('q$\infty$', xy = (V[40],inf[40]), xytext = (V[40],inf[40]))
            ax.set_xlabel('Voltage (mV)')
            show()

# Define subunit kinetics
def alpham(V):
    out = 0.32*(V+37)/(1-exp(-(V+37)*1.0/4))
    return out

def betam(V):
    out = -0.28*(V+10)/(1-exp((V+10)*1.0/5))
    return out

def alphah(V):
    out = 0.128*exp(-(V+33)*1.0/18)
    return out

def betah(V):
    out = 4.0/(1+exp(-(V+10)*1.0/5))
    return out

def alphan(V):
    out = 0.032*(V+35)/(1-exp(-(V+35)*1.0/5))
    return out

def betan(V):
    out = 0.5*exp(-(V+40)*1.0/40)
    return out

def alphap(V):
    out = 1.0/(1+exp(-(V+27.1)/7.2))
    return out

def betap(V):
    out = 21.7-21.0/(3.1+exp(-(V+68.1)*1.0/20.5))
    return out

def alphaq(V):
    out = 1.0/(1+exp((V+32.1)*1.0/5.5))
    return out

def betaq(V):
    out = 105-89.8/(1+exp(-(V+55)*1.0/16.9))
    return out
# End


#Define electrophysiological parameters
mstep=0.01
C = 1
gNa = 50
gK = 25
gCa=0.02
gL = 0.05
ENa = 50
EK = -90
EL = -70
ECa=120



##刺激时间
t1=100 #-ms
t2=320 #-ms320


for Iinj in arange(1.2,7.02,0.582):      #Input current -uA 1.2
    hh = Neuron()
    print  Iinj
    mvoltageArray, mdt,mSpkNum=hh()
    print mSpkNum
    np.savetxt('v.txt',mvoltageArray)
