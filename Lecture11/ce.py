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
        self.s_t=1
    def computer_Isyn(self,t,pre_V,V):
        Vth=-10
        s_inf=tanh((pre_V-Vth)*1.0/10)
        temp=self.s_t
        self.s_t=self.s_t+mstep*((s_inf-self.s_t)*1.0/((1-s_inf)*self.T*self.s_t))
        return self.g_syn*temp*(V-self.E_syn)



class Synaptic_NMDA(object):
    def __init__(self,g_syn,Ts,Tf, Is,If, E_syn):
        self.g_syn=g_syn
        self.Ts=Ts
        self.Tf=Tf
        self.Is=Is
        self.If=If
        self.E_syn=E_syn
        self.preV=[]
        self.msyn=[]
    def computer_Isyn(self,t,t_s,V):
        #self.g_syn*(t-t_s)*1.0/self.T*exp(-(t-t_s)*1.0/self.T)
        #print self.g_syn*(t-t_s)*1.0/self.T*exp(-(t-t_s)*1.0/self.T)
        return self.g_syn*(1.0/(1+0.28*exp(-0.062*V)))*(self.If*exp((t_s-t)*1.0/self.Tf) + self.Is*exp((t_s-t)*1.0/self.Ts))*(V-self.E_syn)
    def func_yita(self, s_Ca):
        yita_0 = 500  #ms
        return yita_0*1.0/(pow(s_Ca,3)+0.0001)
    def computer_Ca(self, Inmda, s_Ca, tao_ca, W):
        gradW = (self.func_oum(s_Ca) - W )/ self.func_yita(s_Ca)
        W  =  W + mstep * gradW
        ##print gradW
        return s_Ca+mstep*(Inmda-s_Ca*1.0/tao_ca),W, self.func_oum(s_Ca),1.0/ self.func_yita(s_Ca)

    def func_sig(self,x,belta):
        return np.exp(belta*x)/(1+np.exp(belta*x))

    def func_oum(self,s_Ca):
        alpha1 = 0.35
        alpha2 = 0.55
        beta1 = 80
        beta2 = 80
        return 0.25 +self.func_sig(s_Ca-alpha2, beta2)-0.25 * self.func_sig(s_Ca-alpha1, beta1)

    def cell_syn(self,t,pre_V,V):
        self.preV.append(pre_V)
        len_=len(self.preV)
        syn_all=0
        if len_>2 and self.preV[-2]>self.preV[-1] and self.preV[-2]>self.preV[-3] and self.preV[-2]>0:
            self.msyn.append(t-mstep)
        if(len(self.msyn)>0):
            for i in range(len(self.msyn)-1,-1,-1):
                I_sy=self.computer_Isyn(t,self.msyn[0],V)
                #print I_sy
            syn_all=syn_all+I_sy
            #print syn_all
        else:
            syn_all=0
        #print syn_all
        return syn_all


class Synaptic_BPAP(object):
    def __init__(self,g_syn,Ts,Tf, Is,If, E_syn):
        self.g_syn=g_syn
        self.Ts=Ts
        self.Tf=Tf
        self.Is=Is
        self.If=If
        self.E_syn=E_syn
        self.preV=[]
        self.msyn=[]
    def computer_Isyn(self,t,t_s,V):
        #self.g_syn*(t-t_s)*1.0/self.T*exp(-(t-t_s)*1.0/self.T)
        #print self.g_syn*(t-t_s)*1.0/self.T*exp(-(t-t_s)*1.0/self.T)
        return self.g_syn*(1.0/(1+0.28*exp(-0.062*V)))*(self.If*exp((t_s-t)*1.0/self.Tf) + self.Is*exp((t_s-t)*1.0/self.Ts))*(V-self.E_syn)
    def func_yita(self, s_Ca):
        yita_0 = 500  #ms
        return yita_0*1.0/(pow(s_Ca,3)+0.0001)
    def computer_Ca(self, Inmda, s_Ca, tao_ca, W):
        gradW = (self.func_oum(s_Ca) - W )/ self.func_yita(s_Ca)
        W  =  W + mstep * gradW
        ##print gradW
        return s_Ca+mstep*(Inmda-s_Ca*1.0/tao_ca),W

    def func_sig(self,x,belta):
        return np.exp(belta*x)/(1+np.exp(belta*x))

    def func_oum(self,s_Ca):
        alpha1 = 0.35
        alpha2 = 0.55
        beta1 = 80
        beta2 = 80
        return 0.25 +self.func_sig(s_Ca-alpha2, beta2)-0.25 * self.func_sig(s_Ca-alpha1, beta1)

    def cell_syn(self,t,pre_V,V):
        self.preV.append(pre_V)
        len_=len(self.preV)
        syn_all=0
        if len_>2 and self.preV[-2]>self.preV[-1] and self.preV[-2]>self.preV[-3] and self.preV[-2]>0:
            self.msyn.append(t-mstep)
        if(len(self.msyn)>0):
            for i in range(len(self.msyn)-1,-1,-1):
                I_sy=self.computer_Isyn(t,self.msyn[0],V)
                #print I_sy
            syn_all=syn_all+I_sy
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
        return (-gL*(self.V-EL)+I)

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
t2=400 #-ms
gmax_ampa=1
gmax_gaba=1
T_ampa=5
E_ampa=0
T_gaba=5
E_gaba=-80

Synaptic_E1=Synaptic(2.56,10,0)
Synaptic_E3=Synaptic(2.56,10,0)
Synaptic_E2=Synaptic(2.56,10,0)
Synaptic_I1=Synaptic(2.56,10,-80)
Synaptic_E4=Synaptic(gmax_ampa*8.0,T_ampa,E_ampa)

Synaptic_nmda1=Synaptic_NMDA(g_syn = -1.0/500,Ts=200,Tf=50, Is=-0.1, If=1.1,E_syn = 130)##g_syn,Ts,Tf, Is,If, E_syn, GNMDA = -1.0/500 ###[uM/(ms*mV)]ENMDA = 130  ##mV
Synaptic_bpap1=Synaptic_BPAP(g_syn = 100.0 ,Ts=25,Tf=3, Is=0.25, If=0.75,E_syn = 100)##g_syn,Ts,Tf, Is,If, E_syn, GNMDA = -1.0/500 ###[uM/(ms*mV)]ENMDA = 130  ##mV
Synaptic_bpap2=Synaptic_BPAP(g_syn = 100.0 ,Ts=25,Tf=3, Is=0.25, If=0.75,E_syn = 100)##g_syn,Ts,Tf, Is,If, E_syn, GNMDA = -1.0/500 ###[uM/(ms*mV)]ENMDA = 130  ##mV
Neuron_A = Neuron()
Neuron_B = Neuron()
Neuron_C = Neuron()
Neuron_D = Neuron()
Neuron_E = Neuron()
A_v=B_v=C_v=D_v=E_v=-65 ##A_skin_Sensory_Neuron, B_tail_Sensory_Neuron, C_interneuron, D_Modulatory_interneuron, E_Motor_neuron
C=1
mtime=arange(t1,t2,0.01)
ter=[]
mCa1=[]
mCa2=[]
s_Ca1=s_Ca2=0
W1=W2=0.25
list_w1=[]
list_w2=[]
list_learn=[]
I_nmda=[]
list_s1=[]
list_s2=[]
IinjA=0
def getI(m_t):
    return -40+40*((1-exp(-t*1.0/50))-1.0/(1+exp(-(t-3800)*1.0/450))+(1-exp(-t*1.0/800)))
for t in mtime:      #Input current -uA
    if t>100 and t<200:
        IinjA=getI(t)
        print IinjA
    elif t>300 and t<400:
        IinjA=0


    A_kv=Neuron_A.solveHHModel(IinjA,A_v)
    syn_AB=Synaptic_I1.computer_Isyn(t,A_v,B_v)
    syn_AC=Synaptic_E1.computer_Isyn(t,A_v,C_v)
    A_v = A_v+mstep*(A_kv)/C
    #Neuron_A.solveHHModel_B()

    C_kv=Neuron_C.solveHHModel(IinjA,C_v)
    C_v = C_v+mstep*(C_kv-syn_AC)/C
    syn_CD=Synaptic_E2.computer_Isyn(t,C_v,D_v)
    #Neuron_C.solveHHModel_B()

    B_kv=Neuron_B.solveHHModel(IinjA,B_v)
    B_v = B_v+mstep*(B_kv-syn_AB)/C
    syn_BE=Synaptic_E3.computer_Isyn(t,B_v,E_v)
    #Neuron_B.solveHHModel_B()

    D_kv=Neuron_D.solveHHModel(IinjA,D_v)
    D_v = D_v+mstep*(D_kv-syn_CD)/C
    #Neuron_D.solveHHModel_B()

    E_kv=Neuron_E.solveHHModel(IinjA,E_v)
    E_v = E_v+mstep*(E_kv-syn_BE)/C
    #Neuron_E.solveHHModel_B()






#plt.plot(list(mtime),ter)
#print Neuron_B.msyn
#print Neuron_B.msyn
#print Neuron_C.msyn
#print len(list(mtime)),len(Neuron_A.mv)
fig = plt.figure()

ax = fig.add_subplot(3,2,1)
ax.set_xlabel('Time (ms)')
ax.set_ylabel('Voltage (mV)')
plt.plot(list(mtime),Neuron_A.mv,'r')
#ax.set_title('figure_1')
ax = fig.add_subplot(3,2,3)
ax.set_xlabel('Time (ms)')
ax.set_ylabel('Voltage (mV)')
plt.plot(list(mtime),Neuron_B.mv,'b')



ax = fig.add_subplot(3,2,4)
ax.set_xlabel('Time (ms)')
ax.set_ylabel('Voltage (mV)')
plt.plot(list(mtime),Neuron_C.mv,'r')

ax = fig.add_subplot(3,2,6)
ax.set_xlabel('Time (ms)')
ax.set_ylabel('SynapsesAB')
plt.plot(list(mtime),Neuron_D.mv,'c')


ax = fig.add_subplot(3,2,5)
ax.set_xlabel('Time (ms)')
ax.set_ylabel('SynapsesCB')
plt.plot(list(mtime),Neuron_E.mv)



show()
#print Neuron_C.mv