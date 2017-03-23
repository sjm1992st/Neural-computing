# coding=utf-8
import numpy as np
import matplotlib.pyplot as plt

class Neuron_IZhi():

    def __init__(self):
        self.v = -65
        self.u  = 0
        self.funcU= 0
        self.vd = -65

    def neu_Mitral(self,I,step):

        if self.v>= -48:
            self.funcU = 20*(self.v+48)

        if self.v< - 48:
            self.funcU = 0

        grad_v = ((self.v+55)*(self.v+50) +0.5*(self.vd - self.v) -self.u +I)/40
        grad_u = 0.4*(self.funcU-self.u)
        grad_vd = 0.0125*(self.v - self.vd)

        self.v  = self.v + step*grad_v
        self.u = self.u + step*grad_u
        self.vd = self.vd +step* grad_vd

        if self.v>=35:
            self.v = -50
            #self.vd = self.vd +step* grad_vd
        #print self.u


        return self.v

    def neu_thalamic(self,I,step):

        grad_v = (0.5*(self.v+60)*(self.v+50) - self.u +I)/20
        grad_u  = 0.05*(7*(self.v +60)-self.u)

        self.v = self.v + step*grad_v
        self.u = self.u +step*grad_u

        if self.v >= 20 - 0.08*self.u:
            self.v = -65 +0.08*self.u
            self.u  = min(self.u+50,530)

        return self.v

    def neu_LS(self,I,step):



        self.grad_v = (0.3*(self.v + 66)*(self.v + 40)+1.2*(self.vd - self.v)-self.u +I)/20
        self.grad_u = 0.17*(5*(self.v+66)-self.u)
        self.grad_vd = 0.01*(self.v -  self.vd)

        if self.v ==30:
            self.v = -45
            self.u = self.u + 100
        else:
            self.v = self.v + step*self.grad_v
            self.u = self.u + step*self.grad_u

        self.vd = self.vd + step* self.grad_vd

        if self.v> 30:
            self.v = 30

        return self.v

class Neuron_Adexp():

    def __init__(self,tao_m,a,tao_w,b,u_r):

        self.v = -60
        #self.u  = -0.01
        self.w = 0
        #self.funcU= 0
        self.tao_m = tao_m
        self.a = a
        self.tao_w = tao_w
        self.b = b
        self.u_r = u_r
        #self.functf = 0

        self.list_v = []
        self.list_t = []
        #self.list_tf = []

    def AdExp(self,I,step):


        v_rest = -70
        R = 0.5 #$m欧姆
        v_th = -50
        delta_ = 2#what i define
        v_peak = 30

        #self.list_v.append(self.v)
        self.grad_v = (- (self.v - v_rest)+delta_ *np.exp((self.v- v_th)/delta_) - R*self.w +R*I)/self.tao_m

        if self.v == v_peak:

            self.grad_w = (self.a*(self.v - v_rest) - self.w +self.b*self.tao_w)/self.tao_w
            self.v = self.u_r
        else:

            self.grad_w = (self.a*(self.v -v_rest) - self.w )/self.tao_w
            self.v = self.v+step*self.grad_v
        '''
        if len(self.list_v)>2 and self.list_v[-2]>self.list_v[-1] and self.list_v[-2]>self.list_v[-3]:

            self.v = self.u_r
            self.grad_w = (self.a*(self.v - v_rest) - self.w +self.b*self.tao_w)/self.tao_w
        '''


        #print self.b*self.tao_w*self.functf
        self.w = self.w + step*self.grad_w

        #if self.w> 20:
            #self.w = self.w +step*self.grad_w

        if  self.v > v_peak:
            self.v  = v_peak

        '''
        if self.v< - 70:
            self.v = -65
        '''
        return self.v

class I_extral():
    def __init__(self,start,end,Iext):
        self.start=start
        self.end=end
        self.Iext=Iext

    def setI(self,t):
        if t<self.start or t>self.end:
            return 0
        else:
            return self.Iext

if __name__ == "__main__":

    #I = 10
    step = 1
    t_start = 0
    t_end = 500 #ms
    step_num =int((t_end - t_start) / step)
    t = np.linspace(t_start, t_end, step_num)
    #v = -65

    Neuron_A = Neuron_IZhi()
    Neuron_B = Neuron_IZhi()
    Neuron_C = Neuron_IZhi()

    Neuron_Bursting = Neuron_Adexp(5.0,-0.5,100.0,7.0,-46.0)
    #Neuron_Bursting = Neuron_Adexp(20,0,30,60,-55)

    Neuron_Irregular = Neuron_Adexp(9.9,-0.5,100,7.0,-46.0)
    Neuron_Transient = Neuron_Adexp(10,1.0,100.0,10.0,-60.0)

    I_A = I_extral(100,400,15)
    I_B = I_extral(100,400,150)
    I_C = I_extral(100,400,9.5)

    I_Bursting = I_extral(100,400,65)
    I_Irregular = I_extral(100,400,65)
    I_Transient = I_extral(100,400,65)

    list_v_A = []
    list_v_B = []
    list_v_C = []

    list_v_Bursting = []
    list_v_Irregular = []
    list_v_Transient = []

    list_i_A = []
    list_i_B = []
    list_i_C = []
    list_i_Bursting = []
    list_i_Irregular = []
    list_i_Transient =[]

    #list_G = []
    for i in t:

        IA = I_A.setI(i)
        list_i_A.append(IA)
        VA = Neuron_A.neu_Mitral(IA,step)

        list_v_A.append(VA)

        IB = I_B.setI(i)
        list_i_B.append(IB)
        VB = Neuron_B.neu_LS(IB,step)
        list_v_B.append(VB)

        IC = I_C.setI(i)
        list_i_C.append(IC)
        VC = Neuron_C.neu_Mitral(IC,step)
        list_v_C.append(VC)

        ID = I_Bursting.setI(i)
        list_i_Bursting.append(ID)
        VD = Neuron_Bursting.AdExp(ID,step)
        print Neuron_Bursting.w
        list_v_Bursting.append(VD)

        IE = I_Irregular.setI(i)
        list_i_Irregular.append(IE)
        VE = Neuron_Irregular.AdExp(IE,step)
        list_v_Irregular.append(VE)

        IF = I_Transient.setI(i)
        list_i_Transient.append(IF)
        VF = Neuron_Transient.AdExp(IF,step)
        list_v_Transient.append(VF)

    fig=plt.figure(1)
    plt.subplot(421)
    plt.plot(t,list_v_A,'r')
    plt.title("neuronA")
    plt.xlabel("t ms")
    plt.ylabel("V mv")
    plt.subplot(423)
    plt.plot(t,list_v_B,'g')
    plt.title("neuronB")
    plt.xlabel("t ms")
    plt.ylabel("V mv")
    plt.subplot(425)
    plt.plot(t,list_v_C,'b')
    plt.title("neuronC")
    plt.xlabel("t ms")
    plt.ylabel("V mv")

    #plt.subplot(4,1,4)
    #plt.plot(t,list_G,'b')
    plt.subplot(427)
    plt.plot(t,list_i_A,'r')
    plt.plot(t,list_i_B,'g')
    plt.plot(t,list_i_C,'b')

    plt.xlabel("t ms")
    plt.ylabel("I pA")

    plt.subplot(422)
    plt.plot(t,list_v_Bursting,'r')
    plt.title("neuron_Bursting")
    plt.xlabel("t ms")
    plt.ylabel("V mv")

    plt.subplot(424)
    plt.plot(t,list_v_Irregular,'g')
    plt.title("neuron_Irregular")
    plt.xlabel("t ms")
    plt.ylabel("V mv")

    plt.subplot(426)
    plt.plot(t,list_v_Transient,'b')
    plt.title("neuron_Transient")
    plt.xlabel("t ms")
    plt.ylabel("V mv")

    plt.subplot(428)
    plt.plot(t,list_i_Bursting,'r')
    plt.plot(t,list_i_Irregular,'g')
    plt.plot(t,list_i_Transient,'b')

    plt.show()