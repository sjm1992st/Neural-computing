# -*- coding: utf-8 -*-
from __future__ import division
from numpy import *
from pylab import *
import matplotlib.pylab as plt
import copy
class Neuron(object):
    def __init__(self):
        self.add_v=[]
        self.cout=0
        self.mv=[]
        self.bw=0
        self.w=0
    def Euler_v(self,I,V):
        #print (gL*(-self.V+EL)+I)
        return (-(V-V_r)+theta*np.exp((V-V_t)*1.0/theta)-R*self.w+R*I)/tao_m

    def Euler_bw(self,V):
        #print (gL*(-self.V+EL)+I)
        return (a*(V-V_r)-self.w+b*tao_w*1*1.0)*1.0/tao_w

    def Euler_w(self,V):
        #print (gL*(-self.V+EL)+I)
        return (a*(V-V_r)-self.w)*1.0/tao_w
    def solveHHModel(self,I,V):
        self.dt=mstep
        self.mv.append(V)
        self.kv1=self.Euler_v(I,V)
        if V==E_peak:
            mw=self.Euler_bw(V)
            V=u_r
        else:

            mw=self.Euler_w(V)
            V = V+mstep*(self.kv1)

        self.w=self.w+mstep*mw
        if V>E_peak:
            V=E_peak
        self.cout=self.cout+1
        return V
        #print kv1,km1



#Define electrophysiological parameters
mstep=1

##刺激时间
t1=0 #-ms
t2=500 #-ms
E_peak=30

A_dv=[]


Neuron_A = Neuron()

V=-60
C=1
mtime=arange(t1,t2,0.01)
ter=[]

tao_w=100
theta=2
R=2
a=-0.5
b=7
V_r=-70
u_r=-46
V_t=-30
tao_m=5
A_flag=0

for t in mtime:      #Input current -uA
    Iinj=2
    print Neuron_A.w

    V=Neuron_A.solveHHModel(Iinj,V)
    #print Neuron_A.V

    #A_dv.append(A_kv)




"""
    if A_v>=V_t:
        A_v=E_peak
        Neuron_A.bw=1
"""


print len(list(mtime)),len(Neuron_A.mv)
fig = plt.figure()
ax = fig.add_subplot(111)
ax.set_xlabel('Time (ms)')
ax.set_ylabel('Voltage (mV)')
plt.plot(list(mtime),Neuron_A.mv)
#ax.set_title('figure_1')

show()

