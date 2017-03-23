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
        self.bw=0
        self.w=0
    def Euler_v(self,I):
        #print (gL*(-self.V+EL)+I)
        return (-(self.V-V_r)+theta*np.exp((self.V-V_t)*1.0/theta)-R*self.w+R*I)

    def Euler_w(self):
        #print (gL*(-self.V+EL)+I)
        return (a*(self.V-V_r)-self.w+b*tao_w*self.bw)*1.0/tao_w

    def solveHHModel(self,I,V):
        self.dt=mstep
        self.V=V

        self.mv.append(self.V)
        self.kv1=self.Euler_v(I)
        self.w=self.w+mstep*self.Euler_w()


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
gmax_ampa=1
gmax_gaba=0.7
T_ampa=5
E_ampa=0
T_gaba=5
E_gaba=-80
E_peak=20

A_dv=[]
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
A_v=B_v=C_v=D_v=E_v=F_v=G_v=-70
C=1
mtime=arange(t1,t2,0.01)
ter=[]

tao_w=30
theta=2
R=20
a=0
b=60
V_r=-55
V_t=-42
tao_m=20
A_flag=0
B_flag=0
C_flag=0
D_flag=0
E_flag=0
F_flag=0
G_flag=0
for t in mtime:      #Input current -uA
    Iinj=6
    print Neuron_A.w
    A_kv=Neuron_A.solveHHModel(Iinj,A_v)
    A_dv.append(A_kv)
    syn_1=Synaptic_1.cell_syn(t,A_v,D_v)
    #print syn_1
    syn_10=Synaptic_10.cell_syn(t,A_v,E_v)
    A_v = A_v+mstep*(A_kv)/tao_m
    if A_flag==1:
        Neuron_A.bw=1
    if A_flag==0:
        Neuron_A.bw=0
    if A_flag==1:
        A_v=V_r
        A_flag=0
    if A_v>V_t:
        A_v=E_peak
        A_flag=1

    #print A_v

    B_kv=Neuron_B.solveHHModel(Iinj,B_v)
    syn_2=Synaptic_2.cell_syn(t,B_v,D_v)
    syn_3=Synaptic_3.cell_syn(t,B_v,E_v)
    syn_4=Synaptic_4.cell_syn(t,B_v,F_v)
    #print 'a',B_v,B_kv,syn_B
    #print t,syn_B
    B_v = B_v+mstep*(B_kv)/tao_m
    if B_flag==1:
        Neuron_B.bw=1
    if B_flag==0:
        Neuron_B.bw=0
    if B_flag==1:
        B_v=V_r
        B_flag=0
    if B_v>V_t:
        B_v=E_peak
        B_flag=1

    #print 'b',B_v


    #print 'B',Neuron_B.mv[-1]
    C_kv=Neuron_C.solveHHModel(Iinj,C_v)
    syn_5=Synaptic_5.cell_syn(t,C_v,E_v)
    syn_6=Synaptic_6.cell_syn(t,C_v,F_v)
    C_v = C_v+mstep*(C_kv)/tao_m
    if C_flag==1:
        Neuron_C.bw=1
    if C_flag==0:
        Neuron_C.bw=0
    if C_flag==1:
        C_v=V_r
        C_flag=0
    if C_v>V_t:
        C_v=E_peak
        C_flag=1


    D_kv=Neuron_D.solveHHModel(0,D_v)
    syn_7=Synaptic_7.cell_syn(t,D_v,G_v)
    D_v = D_v+mstep*(D_kv+syn_1+syn_2)/tao_m
    #print syn_1,syn_2
    if D_flag==1:
        Neuron_D.bw=1
    if D_flag==0:
        Neuron_D.bw=0
    if D_flag==1:
        D_v=V_r
        D_flag=0
    if D_v>V_t:
        D_v=E_peak
        D_flag=1


    E_kv=Neuron_E.solveHHModel(0,E_v)
    syn_8=Synaptic_8.cell_syn(t,E_v,G_v)
    E_v = E_v+mstep*(E_kv+syn_3+syn_5+syn_10)/tao_m
    if E_flag==1:
        Neuron_E.bw=1
    if E_flag==0:
        Neuron_E.bw=0
    if E_flag==1:
        E_v=V_r
        E_flag=0
    if E_v>V_t:
        E_v=E_peak
        E_flag=1


    F_kv=Neuron_F.solveHHModel(0,F_v)
    syn_9=Synaptic_9.cell_syn(t,F_v,G_v)
    F_v = F_v+mstep*(F_kv+syn_4+syn_6)/tao_m
    if F_flag==1:
        Neuron_F.bw=1
    if F_flag==0:
        Neuron_F.bw=0
    if F_flag==1:
        F_v=V_r
        F_flag=0
    if F_v>V_t:
        F_v=E_peak
        F_flag=1


    G_kv=Neuron_G.solveHHModel(0,G_v)
    G_v = G_v+mstep*(G_kv+syn_7+syn_8+syn_9)/tao_m
    if G_flag==1:
        Neuron_G.bw=1
    if G_flag==0:
        Neuron_G.bw=0
    if G_flag==1:
        G_v=V_r
        G_flag=0
    if G_v>V_t:
        G_v=E_peak
        G_flag=1


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

