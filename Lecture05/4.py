# -*- coding: utf-8 -*-
from __future__ import division
from numpy import *
from pylab import *
import matplotlib.pylab as plt
import copy
class Synaptic(object):
    def __init__(self,g_syn,T,E_syn):
        self.g_syn=g_syn
        self.T=T
        self.E_syn=E_syn
        self.preV=[]
        self.msyn=[]
    def computer_Isyn(self,t,t_s,V):
        #self.g_syn*(t-t_s)*1.0/self.T*exp(-(t-t_s)*1.0/self.T)
        #print self.g_syn*(t-t_s)*1.0/self.T*exp(-(t-t_s)*1.0/self.T)
        return self.g_syn*(t-t_s)*1.0/self.T*exp(-(t-t_s)*1.0/self.T)*(V-self.E_syn)
    def cell_syn(self,t,pre_V,V):
        self.preV.append(pre_V)
        len_=len(self.preV)
        syn_all=0
        if len_>2 and self.preV[-2]>self.preV[-1] and self.preV[-2]>self.preV[-3] and self.preV[-2]>0:
            self.msyn.append(t-mstep)
        if(len(self.msyn)>0):
            for i in range(len(self.msyn)-1,-1,-1):
                I_sy=self.computer_Isyn(t,self.msyn[i],V)
                #print I_sy
                if abs(I_sy)<0.01:
                    break
                syn_all=syn_all-I_sy
            #print syn_all
        else:
            syn_all=0
        #print syn_all
        return syn_all
class Neuron(object):
    def __init__(self):
        self.add_v=[]
        self.cout=0
        self.mv=[]
    def Euler_parameter(self,p_va,p_vb,kv):
         return p_va*(1-kv) - p_vb*kv
    def Euler_v(self,I):
        return (-gNa*self.m*self.m*self.m*self.h*(self.V-ENa)-gK*self.n*self.n*self.n*self.n*(self.V-EK)-gL*(self.V-EL)+I)

    def solveHHModel(self,I,V):
        self.dt=mstep
        if self.cout==0:
            self.V = -65
            self.m = 0.053
            self.h = 0.596
            self.n = 0.318
        else:
            self.V=V
        self.mv.append(self.V)
        self.kv1=self.Euler_v(I)
        self.km1=self.dt*self.Euler_parameter(alpham(self.V),betam(self.V),self.m)
        self.kh1=self.dt*self.Euler_parameter(alphah(self.V),betah(self.V),self.h)
        self.kn1=self.dt*self.Euler_parameter(alphan(self.V),betan(self.V),self.n)
        self.cout=self.cout+1
        return self.kv1
        #print kv1,km1
    def solveHHModel_B(self):
        self.m = self.m+self.km1
        self.h = self.h+self.kh1
        self.n = self.n+self.kn1
        #print V



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
t1=0 #-ms
t2=150 #-ms
gmax_ampa=1
gmax_gaba=1
T_ampa=5
E_ampa=0
T_gaba=5
E_gaba=-80

Synaptic_E1=Synaptic(gmax_ampa,T_ampa,E_ampa)
Synaptic_E2=Synaptic(gmax_ampa*1.0/3,T_ampa,E_ampa)
Synaptic_I1=Synaptic(gmax_gaba,T_gaba,E_gaba)
Synaptic_I2=Synaptic(gmax_gaba*1.0/2,T_gaba,E_gaba)
Neuron_A = Neuron()
Neuron_B = Neuron()
Neuron_C = Neuron()
A_v=B_v=C_v=-65
C=1
mtime=arange(t1,t2,0.01)
ter=[]
for t in mtime:      #Input current -uA
    IinjA=15
    A_kv=Neuron_A.solveHHModel(IinjA,A_v)
    syn_B=Synaptic_I2.cell_syn(t,A_v,C_v)
    A_v = A_v+mstep*(A_kv)/C
    Neuron_A.solveHHModel_B()
    #print 'A',Neuron_A.mv[-1]
    ter.append(syn_B)
    IinjB=10
    B_kv=Neuron_B.solveHHModel(IinjB,B_v)
    syn_C=Synaptic_E2.cell_syn(t,B_v,C_v)
    #print 'a',B_v,B_kv,syn_B
    #print t,syn_B
    B_v = B_v+mstep*(B_kv)/C
    #print 'b',B_v
    Neuron_B.solveHHModel_B()

    #print 'B',Neuron_B.mv[-1]
    C_kv=Neuron_C.solveHHModel(0,C_v)
    C_v = C_v+mstep*(C_kv+syn_B+syn_C)/C
    Neuron_C.solveHHModel_B()
   # print 'C',Neuron_C.mv[-1]
#plt.plot(list(mtime),ter)
#print Neuron_B.msyn
#print Neuron_B.msyn
#print Neuron_C.msyn
print len(list(mtime)),len(Neuron_A.mv)
fig = plt.figure()
ax = fig.add_subplot(311)
ax.set_xlabel('Time (ms)')
ax.set_ylabel('Voltage (mV)')
plt.plot(list(mtime),Neuron_A.mv)
#ax.set_title('figure_1')
ax = fig.add_subplot(312)
ax.set_xlabel('Time (ms)')
ax.set_ylabel('Voltage (mV)')
plt.plot(list(mtime),Neuron_B.mv)

ax = fig.add_subplot(313)
ax.set_xlabel('Time (ms)')
ax.set_ylabel('Voltage (mV)')
plt.plot(list(mtime),Neuron_C.mv)
show()
#print Neuron_C.mv