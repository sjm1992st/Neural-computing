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
    def Euler_v(self,I):
        #print (gL*(-self.V+EL)+I)
        return (-(self.V-V_r)+theta*np.exp((self.V-V_t)*1.0/theta)-R*self.w+R*I)

    def Euler_w(self):
        #print (gL*(-self.V+EL)+I)
        return (a*(self.V-V_r)-self.w+b*tao_w*self.bw*1.0)*1.0/tao_w

    def solveHHModel(self,I,V):
        self.dt=mstep
        self.V=V

        self.mv.append(self.V)

        if self.bw==1:
            self.w=self.w+b
        self.kv1=self.Euler_v(I)
        self.w=self.w+mstep*self.Euler_w()


        self.cout=self.cout+1
        return self.kv1
        #print kv1,km1



#Define electrophysiological parameters
mstep=0.01

##刺激时间
t1=0 #-ms
t2=500 #-ms
E_peak=38

A_dv=[]


Neuron_A = Neuron()

A_v=-70
C=1
mtime=arange(t1,t2,0.01)
ter=[]

tao_w=100
theta=2
R=2
a=1
b=10
V_r=-60
V_t=-38
tao_m=910
A_flag=0

for t in mtime:      #Input current -uA
    Iinj=65
    print Neuron_A.w
    A_kv=Neuron_A.solveHHModel(Iinj,A_v)
    #print Neuron_A.V
    A_dv.append(A_kv)
    A_v = A_v+mstep*(A_kv)/tao_m
    if A_flag==1:
        Neuron_A.bw=1
    if A_flag==0:
        Neuron_A.bw=0
    if A_flag==1:
        A_v=V_r
        A_flag=0
    if A_v>=V_t:
        A_v=E_peak
        A_flag=1


print len(list(mtime)),len(Neuron_A.mv)
fig = plt.figure()
ax = fig.add_subplot(111)
ax.set_xlabel('Time (ms)')
ax.set_ylabel('Voltage (mV)')
plt.plot(list(mtime),Neuron_A.mv)
#ax.set_title('figure_1')

show()

