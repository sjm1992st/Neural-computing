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
        self.u=0
        self.vd=-75
    def Euler_v(self,I):
        #print (gL*(-self.V+EL)+I)
        return (K*(self.V-V_r)*(self.V-V_t)-self.u+I+1.2*(self.vd-self.V))

    def Euler_U(self):
        #print (gL*(-self.V+EL)+I)
        return (a*(b*(self.V-V_r)-self.u))
    def Euler_Vd(self):
        #print (gL*(-self.V+EL)+I)
        return 0.01*(self.V-self.vd)

    def solveHHModel(self,I,V):
        self.dt=mstep
        self.V=V
        self.mv.append(self.V)
        self.kv1=self.Euler_v(I)
        self.vd=self.vd+mstep*self.Euler_Vd()
        self.u=self.u+mstep*self.Euler_U()


        self.cout=self.cout+1
        return self.kv1
        #print kv1,km1



#Define electrophysiological parameters
mstep=0.01
C = 1
gL = 0.3
EL = -40



##刺激时间
t1=0 #-ms
t2=250 #-ms

E_peak=30

A_dv=[]


Neuron_A = Neuron()

A_v=-75
C=20
mtime=arange(t1,t2,0.01)
ter=[]
K=0.3
a=0.17
b=5
c=-45
d=100
V_r=-66
V_t=-40
for t in mtime:      #Input current -uA
    Iinj=150
    A_kv=Neuron_A.solveHHModel(Iinj,A_v)
    A_dv.append(A_kv)
    A_v = A_v+mstep*(A_kv)/C
    if A_v>=E_peak:
        A_v=c
        Neuron_A.u=Neuron_A.u+d
    #print A_v

print len(list(mtime)),len(Neuron_A.mv)
fig = plt.figure()
ax = fig.add_subplot(111)
ax.set_xlabel('Time (ms)')
ax.set_ylabel('Voltage (mV)')
plt.plot(list(mtime),Neuron_A.mv)
#ax.set_title('figure_1')

show()

