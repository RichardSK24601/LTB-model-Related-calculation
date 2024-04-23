# -*- coding: utf-8 -*-
"""
Created on Mon Apr 15 13:27:40 2024

@author: 87588
"""

from threading import Thread
from multiprocessing import Process
import math
import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
from scipy.integrate import quad
from tqdm import trange, tqdm
from scipy.integrate import odeint
from tqdm import trange, tqdm

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
def intg(r,t):
    return rho(r,t)*r**2

def H_frw(t):
    return np.sqrt((4/9)*rho_frw(t))
def rho_bar(r,t):
    temp=quad(intg, 0, r,args=t)[0]
    temp=temp*3/r**3
    return temp

def M(r):
    r=np.abs(r)
    t=1
    temp=rho_bar(r,t)*r**3*4/3*np.pi
    return temp
def dM(r):
    dr=0.001
    return (M(r+dr/2)-M(r-dr/2))/dr

def E(r):
    r=np.abs(r)
    t=1
    temp=1/2*H_frw(t)**2*r**2-1/(6*np.pi)*M(r)/r
    return temp

def dE(r):
    dr=0.001
    return (E(r+dr/2)-E(r-dr/2))/dr

def dYdt1(y,t,r):
    temp=M(r)/(3*np.pi)*(1/y)+2*E(r)
    return np.sqrt(temp)

def Y(r,t):
    r1=np.abs(r)
    ttemp=np.linspace(-0.8,t,200)
    sol = odeint(dYdt1, r1, ttemp, args=(r1,))
    return sol[len(sol)-1][0]

def dYdt(r,t):
    step=0.0001
    return (Y(r,t+step/2)-Y(r,t-step/2))/step

def dYdr(r,t):
    step=0.0001
    return (Y(r+step/2,t)-Y(r-step/2,t))/step

def dYdrdt(r,t):
    step=0.0001
    return (dYdt(r+step/2,t)-dYdt(r-step/2,t))/step

def W(r):
    return np.sqrt(2*E(np.absolute(r))+1)

def fdr(r,t,z,c):
    temp0=(z+1)**2-c**2/Y(r,t)**2
    temp1=(W(r)/dYdr(r,t))*np.sqrt(np.abs(temp0))
    return temp1

def infdr(r,t,z,c):
    temp0=(z+1)**2-c**2/Y(r,t)**2
    temp1=(W(r)/dYdr(r,t))*np.sqrt(np.abs(temp0))
    return -temp1
'''
def fdr(r,t,z,c):
    temp0=(z+1)**2-c**2/Y(r,t)**2
    temp1=(W(r)/dYdr(r,t))*np.sqrt(temp0)
    return temp1
'''

def intersec(l,r,z,t,phi):
    rh=0.042
    phi=np.array(phi)
    r=np.array(r)
    f=r**2+3*rh**2-4*r*rh*np.cos(phi)
    for i in range(len(f)-1):
        if(f[i]>0 and f[i+1]<0):
            return l[i],r[i],z[i],t[i],phi[i]
        if(f[i]<0 and f[i+1]>0):
            return l[i],r[i],z[i],t[i],phi[i]
    

def geo(c1,rint,zint,tint,phiint,rcut):
    steps=[]
    rArray=[]
    tArray=[]
    zArray=[]
    phiArray=[]
    flip1=0
    c3=c1
    def evo(z,t,r,phi,flip,c4):
        #global c
        c=c4
        step=0.0001
        if(flip==0):
            if((z+1)**2-c**2/Y(r,t)**2<0):
                print('flip')
                flip=1
                phi=np.pi-phi
        if (flip==0):
            tempr=r-step*fdr(r,t,z,c)
            tempz=z-step*(-(dYdrdt(r,t)/dYdr(r,t))*((z+1)**2-c**2/Y(r,t)**2)-c**2*(dYdt(r,t)/Y(r,t)**3))
            tempphi=phi-step*c**2/Y(r,t)**2
            tempt=t-step*(z+1)
        if(flip==1):
            tempr=r+step*fdr(r,t,z,c)
            tempz=z-step*(-(dYdrdt(r,t)/dYdr(r,t))*((z+1)**2-c**2/Y(r,t)**2)-c**2*(dYdt(r,t)/Y(r,t)**3))
            tempphi=phi-step*c**2/Y(r,t)**2
            tempt=t-step*(z+1)
        return tempr,tempt,tempz,tempphi,flip
    
    for i in tqdm(range(0,2300)):
        rArray.append(rint)
        tArray.append(tint)
        zArray.append(zint)
        phiArray.append(phiint)
        steps.append(i*0.01)
        rint, tint, zint, phiint, flip1=evo(zint,tint,rint,phiint,flip1,c3)
        if (rint>rcut):
            break
    return steps,rArray,tArray,zArray,phiArray
#%%
#testbench
'''
rint=0.042
zint=0
tint=-0.25
phiint=np.pi
c1=0.01
'''
if __name__ == "__main__":
    rh=0.042
    ri=0.042
    zi=0
    ti=0
    phii=np.pi
    c1=0
    rc=0.043
    l,r,t,z,phi=geo(c1,ri,zi,ti,phii,rc)
    plt.figure()
    plt.plot(r,z)
    plt.figure()
    plt.plot(r,phi)
    plt.figure()
    plt.plot(r,t)
    