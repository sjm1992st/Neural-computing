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
        self.vd=-70
    def Euler_v(self,I):
        #print (gL*(-self.V+EL)+I)
        return (K*(self.V-V_r)*(self.V-V_t)-self.u+I+0.5*(self.vd-self.V))

    def Euler_U(self,UU):
        #print (gL*(-self.V+EL)+I)
        return (a*(UU-self.u))
    def Euler_Vd(self):
        #print (gL*(-self.V+EL)+I)
        return 0.0125*(self.V-self.vd)

    def solveHHModel(self,I,V,UU):
        self.dt=mstep
        self.V=V
        self.mv.append(self.V)
        self.kv1=self.Euler_v(I)
        self.vd=self.vd+mstep*self.Euler_Vd()
        self.u=self.u+mstep*self.Euler_U(UU)


        self.cout=self.cout+1
        return self.kv1
        #print kv1,km1



#Define electrophysiological parameters
mstep=0.01


##刺激时间
t1=0 #-ms
t2=300 #-ms

E_peak=30

A_dv=[]


Neuron_A = Neuron()

A_v=-70
C=40
mtime=arange(t1,t2,0.01)
ter=[]
K=1
a=0.4
b=5
c=-50
d=100
V_r=-55
V_t=-50
for t in mtime:      #Input current -uA
    Iinj=15
    if A_v>=-48:
        UU=20*(A_v+48)
    if A_v< -48:
        UU=0
    if A_v>=35:
        A_v=c
    A_kv=Neuron_A.solveHHModel(Iinj,A_v,UU)
    A_dv.append(A_kv)
    A_v = A_v+mstep*(A_kv)/C
    #print A_v

print len(list(mtime)),len(Neuron_A.mv)
fig = plt.figure()
ax = fig.add_subplot(111)
ax.set_xlabel('Time (ms)')
ax.set_ylabel('Voltage (mV)')
plt.plot(list(mtime),Neuron_A.mv)
#ax.set_title('figure_1')

show()

