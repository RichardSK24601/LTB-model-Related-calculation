# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

import numpy as np
from scipy.integrate import quad
import matplotlib.pyplot as plt
def rho_frw(t):
    return 25
def rho(r,t):
    r_h=0.042
    A=50.59
    r_M=0.037
    sig=r_h/10
    eps=0.0025

    if r> r_h:
        return rho_frw(t)
    else:
        return A*np.exp(-1*(r-r_M)**2/(2*sig**2))+eps
def Y(r,t):
    return r
def intg(r,t):
    return rho(r,t)*r**2

def H_frw(t):
    return np.sqrt((4/9)*rho_frw(t))
def rho_bar(r,t):
    temp=quad(intg, 0, r,args=t)[0]
    temp=temp*3/r**3
    return temp

def M(r):
    t=1
    temp=rho_bar(r,t)*r**3*4/3*np.pi
    return temp
def E(r):
    t=1
    temp=1/2*H_frw(t)**2*r**2-1/(6*np.pi)*M(r)/r
    return temp
#%%
#
r=np.linspace(0, 0.05,100)
def E1(r):
    t=1
    temp=1/2*H_frw(t)**2*r**2
    return temp

def E2(r):
    t=1
    temp=1/(6*np.pi)*M(r)/r
    return temp

e1=[]
e2=[]
for i in range(0,len(r)):
    e1.append(E1(r[i]))
    e2.append(E2(r[i]))

e1=np.array(e1)
e2=np.array(e2)
e3=e1-e2
plt.figure()
plt.plot(r,e1)
plt.plot(r,e2)
plt.plot(r,e3)
#%%
#Fig 2
r=np.linspace(0, 0.05,10000)
t=1
I=[]
for i in range(0,len(r)):
    #temp=quad(intg, 0, r[i],args=t)[0]
    #temp=temp*3/r[i]**3
    #I.append(temp)
    I.append(rho_bar(r[i],t))
I=np.array(I)
rplot=[]
for i in range(0,len(r)):
    rplot.append(rho(r[i], t))

plt.figure()
plt.plot(r,rplot)
plt.plot(r,I)
#%%
#Fig.4
r=np.linspace(0, 0.05,10000)
E1=[]
for i in range (0,len(r)):
    E1.append(E(r[i]))
plt.figure()
plt.plot(r,E1)