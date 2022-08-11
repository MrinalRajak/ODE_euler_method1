

#comparison of different methods of solving ODE's by
#considering L-R ckt as test case
#L(dx/dy)+XR=E
#X-->current,E-->emf,R-->Resistance,L-->Inductance,t-->time

import numpy as np
import matplotlib.pyplot as plt
def f(x,t):
     return (E-x*R)/L
E=1.0
R=1.0
L=1.0

def Euler(f,tt):
    t=np.arange(tt[0],tt[1],tt[2])
    x=np.empty(len(t))
    x[0]=tt[3]
    for k in range(0,len(t)-1):
        k1=f(x[k],t[k])
        x[k+1]=x[k]+tt[2]*k1
    return np.array([t,x])

def Mid_Euler(f,tt):
    t=np.arange(tt[0],tt[1],tt[2])
    x=np.empty(len(t))
    x[0]=tt[3]
    for k in range(0,len(t)-1):
        k1=f(x[k],t[k])
        k2=f(x[k]+tt[2]*k1/2.0,t[k]+tt[2]/2.0)
        x[k+1]=x[k]+tt[2]*k2
    return np.array([t,x])

def Heun(f,tt):
    t=np.arange(tt[0],tt[1],tt[2])
    x=np.empty(len(t))
    x[0]=tt[3]
    for k in range(0,len(t)-1):
        k1=f(x[k],t[k])
        k2=f(x[k]+tt[2]*k1,t[k]+tt[2])
        x[k+1]=x[k]+tt[2]*(k1+k2)/2.0
    return np.array([t,x])

def Rk4(f,tt):
    t=np.arange(tt[0],tt[1],tt[2])
    x=np.empty(len(t))
    x[0]=tt[3]
    for k in range(0,len(t)-1):
        k1=tt[2]*f(x[k],t[k])
        k2=tt[2]*f(x[k]+tt[2]*k1/2.0,t[k]+tt[2]/2.0)
        k3=tt[2]*f(x[k]+tt[2]*k2/2.0,t[k]+tt[2]/2.0)
        k4=tt[2]*f(x[k]+tt[2]*k3/2.0,t[k]+tt[2]/2.0)
        x[k+1]=x[k]+(k1+2*k2+2*k3+k4)/6.0
    return np.array([t,x])

def exact(f,tt):
    t=np.arange(tt[0],tt[1],tt[2])
    exact=R*(1-np.exp((-R*t)/L))
    return np.array([t,exact])

tt=[0.0,10.0,0.01,0.0] # init time,fin time,stepsize,init
t,Xe=Euler(f,tt)
t,Xme=Mid_Euler(f,tt)
t,Xhe=Heun(f,tt)
t,Xrk4=Rk4(f,tt)
t,exact=exact(f,tt)
error_Euler=np.abs(exact-Xe)
error_Mid=np.abs(exact-Xme)
error_heun=np.abs(exact-Xhe)
error_Rk4=np.abs(exact-Xrk4)
fig,axrr=plt.subplots(4,2)
axrr[0,0].plot(t,Xe)
axrr[0,0].set_title(r'$Euler$',fontsize=15)
axrr[0,1].plot(t,error_Euler)
axrr[0,1].set_title(r'$ErrorEuler$',fontsize=15)
axrr[1,0].plot(t,Xme)
axrr[1,0].set_title(r'$MidEuler$',fontsize=15)
axrr[1,1].plot(t,error_Mid)
axrr[1,1].set_title(r'$ErrorMidEuler$',fontsize=15)
axrr[2,0].plot(t,Xhe)
axrr[2,0].set_title(r'$HeunEuler$',fontsize=15)
axrr[2,1].plot(t,error_heun)
axrr[2,1].set_title(r'$ErrorHeunEuler$',fontsize=15)
axrr[3,0].plot(t,Xrk4)
axrr[3,0].set_title(r'$Rk4Euler$',fontsize=15)
axrr[3,1].plot(t,error_Rk4)
axrr[3,1].set_title(r'$ErrorRk4$',fontsize=15)
fig.subplots_adjust(hspace=0.5,wspace=0.3)
plt.show()
        
        





































































    
