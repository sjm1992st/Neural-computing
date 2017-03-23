# -*- coding: utf-8 -*-
from __future__ import division
from numpy import *
from pylab import *
import matplotlib.pylab as plt

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


hh = Neuron()
mvoltageArray, mdt,mSpkNum=hh()
print mSpkNum

############np.linspace(-20,80,100)
v_1=np.arange(-76,15,0.1)
n_1=np.arange(0,1,0.05)
n_list0=np.zeros([len(v_1),1])
n_list1=np.zeros([len(v_1),1])
I=4.51
kk=0
for v_e in v_1:
    n_list0[kk,0]=(n_inf(v_e))
    n_=(I-gNa*m_inf(v_e)*(v_e-ENa)-gL*(v_e-EL))*1.0/gK/(v_e-EK)
    n_list1[kk,0]=n_
    kk=kk+1
#print v_1,n_list
nk=np.argmin(abs(n_list0-n_list1))
print nk
plt.plot(v_1[nk],n_list0[nk],'go')
plt.plot(v_1,n_list0,'r')
plt.plot(v_1,n_list1)

####vector
deta=0.1
v_2=np.arange(-80,15,5)
n_2=np.arange(0,1,0.05)
In=I
for ni in range(len(n_2)):
    for vj in range(len(v_2)):
        if ni>0 and ni<len(n_2)-1 and vj>0 and vj<len(v_2)-1:

            dv=(-gNa*m_inf(v_2[vj])*(v_2[vj]-ENa)-gK*n_2[ni]*(v_2[vj]-EK)-gL*(v_2[vj]-EL)+In)/C
            dn=n_inf(v_2[vj])-n_2[ni]
            """
            v_pre=v_2[vj]-deta
            a=dn*1.0/dv
            print a
            n_pre=n_2[ni]-deta*a
            v_last=v_2[vj]+deta
            n_last=n_2[ni]+deta*a
            #print v_pre,n_pre
            #plt.plot([v_pre,v_last],[n_pre,n_last],'k')
            print a
            if a>0:
                plt.annotate('',xy=(v_last,n_last),xytext=(v_pre,n_pre),arrowprops=dict(arrowstyle="->",connectionstyle="arc3"))

           """
            ax=dv*0.01
            ay=dn*0.01
            #print dv,dn
            plt.annotate('',xy=(v_2[vj],n_2[ni]),xytext=(v_2[vj]+ax,n_2[ni]+ay),arrowprops=dict(arrowstyle="<-",connectionstyle="arc3"))
            #plt.plot([v_2[vj],v_2[vj]+ax],[n_2[ni],n_2[ni]+ay],'k')

##limit cycle
xpos=-20
ypos=0.2
kl=0

while(kl<1500):
    dv=(-gNa*m_inf(xpos)*(xpos-ENa)-gK*ypos*(xpos-EK)-gL*(xpos-EL)+In)/C
    dn=n_inf(xpos)-ypos
    xnpos=xpos+dv*0.01
    ynpos=ypos+dn*0.01
    plt.plot([xpos,xnpos],[ypos,ynpos],'g')
    if kl%100==0 and kl>0:
        plt.annotate('',xy=(xpos,ypos),xytext=(xnpos,ynpos),arrowprops=dict(arrowstyle="<-",connectionstyle="arc3"))
    kl=kl+1
    xpos=xnpos
    ypos=ynpos
plt.show()

