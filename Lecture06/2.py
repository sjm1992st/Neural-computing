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
                if abs(I_sy)<0.001:
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
    def Euler_v(self,I):
        #print (gL*(-self.V+EL)+I)
        return (gL*(-self.V+E_peak)**2/st+I)

    def solveHHModel(self,I,V):
        self.dt=mstep
        if self.cout==0:
            self.V = 0
            self.m = 0.053
            self.h = 0.596
            self.n = 0.318
        else:
            self.V=V

        self.mv.append(self.V)
        self.kv1=self.Euler_v(I)

        self.cout=self.cout+1
        return self.kv1
        #print kv1,km1



#Define electrophysiological parameters
mstep=0.01
C = 1
gNa = 120
gK = 36
gL = 0.3
ENa = 50
EK = -77
EL = -40



##刺激时间
t1=0 #-ms
t2=150 #-ms
gmax_ampa=0.5
gmax_gaba=1
T_ampa=5
E_ampa=24
T_gaba=5
E_gaba=-80
E_peak=20
st=6

Synaptic_1=Synaptic(gmax_ampa,T_ampa,E_ampa)
Synaptic_2=Synaptic(gmax_ampa,T_ampa,E_ampa)
Synaptic_3=Synaptic(gmax_ampa,T_ampa,E_ampa)
Synaptic_4=Synaptic(gmax_ampa,T_ampa,E_ampa)
Synaptic_6=Synaptic(gmax_ampa,T_ampa,E_ampa)
Synaptic_7=Synaptic(gmax_ampa,T_ampa,E_ampa)
Synaptic_8=Synaptic(gmax_ampa,T_ampa,E_ampa)
Synaptic_9=Synaptic(gmax_ampa,T_ampa,E_ampa)
Synaptic_5=Synaptic(gmax_gaba*1.0/10,T_gaba,E_gaba)
Synaptic_10=Synaptic(gmax_gaba*1.0/10,T_gaba,E_gaba)

Neuron_A = Neuron()
Neuron_B = Neuron()
Neuron_C = Neuron()
Neuron_D = Neuron()
Neuron_E = Neuron()
Neuron_F = Neuron()
Neuron_G = Neuron()
A_v=B_v=C_v=D_v=E_v=F_v=G_v=-80
C=1
mtime=arange(t1,t2,0.01)
ter=[]
for t in mtime:      #Input current -uA
    Iinj=1
    A_kv=Neuron_A.solveHHModel(Iinj,A_v)
    syn_1=Synaptic_1.cell_syn(t,A_v,D_v)
    #print syn_1
    syn_10=Synaptic_10.cell_syn(t,A_v,E_v)
    A_v = A_v+mstep*(A_kv)/C
    if A_v>E_peak:
        A_v=EL
    #print A_v

    B_kv=Neuron_B.solveHHModel(Iinj,B_v)
    syn_2=Synaptic_2.cell_syn(t,B_v,D_v)
    syn_3=Synaptic_3.cell_syn(t,B_v,E_v)
    syn_4=Synaptic_4.cell_syn(t,B_v,F_v)
    #print 'a',B_v,B_kv,syn_B
    #print t,syn_B
    B_v = B_v+mstep*(B_kv)/C
    if B_v>E_peak:
        B_v=EL
    #print 'b',B_v


    #print 'B',Neuron_B.mv[-1]
    C_kv=Neuron_C.solveHHModel(Iinj,C_v)
    syn_5=Synaptic_5.cell_syn(t,C_v,E_v)
    syn_6=Synaptic_6.cell_syn(t,C_v,F_v)
    C_v = C_v+mstep*(C_kv)/C
    if C_v>E_peak:
        C_v=EL

    D_kv=Neuron_D.solveHHModel(0,D_v)
    syn_7=Synaptic_7.cell_syn(t,D_v,G_v)
    D_v = D_v+mstep*(D_kv+syn_1+syn_2)/C
    #print syn_1,syn_2
    if D_v>E_peak:
        D_v=EL

    E_kv=Neuron_E.solveHHModel(0,E_v)
    syn_8=Synaptic_8.cell_syn(t,E_v,G_v)
    E_v = E_v+mstep*(E_kv+syn_3+syn_5+syn_10)/C
    if E_v>E_peak:
        E_v=EL

    F_kv=Neuron_F.solveHHModel(0,F_v)
    syn_9=Synaptic_9.cell_syn(t,F_v,G_v)
    F_v = F_v+mstep*(F_kv+syn_4+syn_6)/C
    if F_v>E_peak:
        F_v=EL

    G_kv=Neuron_G.solveHHModel(0,G_v)
    G_v = G_v+mstep*(G_kv+syn_7+syn_8+syn_9)/C
    if G_v>E_peak:
        G_v=EL

   # print 'C',Neuron_C.mv[-1]
#plt.plot(list(mtime),ter)
#print Neuron_B.msyn
#print Neuron_B.msyn
#print Neuron_C.msyn
print len(list(mtime)),len(Neuron_A.mv)
fig = plt.figure()
ax = fig.add_subplot(431)
ax.set_xlabel('Time (ms)')
ax.set_ylabel('Voltage (mV)')
plt.plot(list(mtime),Neuron_A.mv)
#ax.set_title('figure_1')

ax = fig.add_subplot(432)
ax.set_xlabel('Time (ms)')
ax.set_ylabel('Voltage (mV)')
plt.plot(list(mtime),Neuron_B.mv)

ax = fig.add_subplot(433)
ax.set_xlabel('Time (ms)')
ax.set_ylabel('Voltage (mV)')
plt.plot(list(mtime),Neuron_C.mv)

ax = fig.add_subplot(434)
ax.set_xlabel('Time (ms)')
ax.set_ylabel('Voltage (mV)')
plt.plot(list(mtime),Neuron_D.mv)

ax = fig.add_subplot(435)
ax.set_xlabel('Time (ms)')
ax.set_ylabel('Voltage (mV)')
plt.plot(list(mtime),Neuron_E.mv)

ax = fig.add_subplot(436)
ax.set_xlabel('Time (ms)')
ax.set_ylabel('Voltage (mV)')
plt.plot(list(mtime),Neuron_F.mv)

ax = fig.add_subplot(437)
ax.set_xlabel('Time (ms)')
ax.set_ylabel('Voltage (mV)')
plt.plot(list(mtime),Neuron_G.mv)
show()
#print Neuron_C.mv
