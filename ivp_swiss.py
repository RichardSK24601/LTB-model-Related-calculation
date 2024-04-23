# -*- coding: utf-8 -*-
"""
Created on Tue Mar 19 16:45:12 2024

@author: 87588
"""

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
#%%
c=0.001
for k in range(0,20):
    steps=[]
    rArray=[]
    tArray=[]
    zArray=[]
    phiArray=[]
    rint=0.042
    zint=0
    tint=0
    phiint=np.pi
    flip=0
    def evo(z,t,r,phi):
        global flip
        global c
        step=0.0001
        if(flip==0):
            if((z+1)**2-c**2/Y(r,t)**2<0):
                print('flip')
                flip=1
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
        return tempr,tempt,tempz,tempphi

    for i in tqdm(range(0,2300)):
        rArray.append(rint)
        tArray.append(tint)
        zArray.append(zint)
        phiArray.append(phiint)
        steps.append(i*0.01)
        rint, tint, zint, phiint=evo(zint,tint,rint,phiint)
        
    string="c=%f"%c
    plt.figure()
    plt.plot(rArray,zArray,label=string) 
    plt.xlabel('r')
    plt.ylabel('z')
    plt.legend()
    plt.savefig('C:/Users/87588/Documents/RG Cosmo/Code/Figure/RvsZ%f.png'%c)
    plt.show()

    plt.figure()
    plt.plot(steps,zArray,label=string) 
    plt.xlabel('l')
    plt.ylabel('z')
    plt.legend()
    plt.savefig('C:/Users/87588/Documents/RG Cosmo/Code/Figure/Zvsl%f.png'%c)
    plt.show()

    plt.figure()
    plt.plot(steps,tArray,label=string) 
    plt.xlabel('l')
    plt.ylabel('t')
    plt.legend()
    plt.savefig('C:/Users/87588/Documents/RG Cosmo/Code/Figure/Tvsl%f.png'%c)
    plt.show()

    plt.figure()
    plt.plot(steps,rArray,label=string) 
    plt.xlabel('l')
    plt.ylabel('r')
    plt.legend()
    plt.savefig('C:/Users/87588/Documents/RG Cosmo/Code/Figure/Rvsl%f.png'%c)
    plt.show()

    plt.figure()
    plt.plot(steps,phiArray,label=string) 
    plt.xlabel('l')
    plt.ylabel('phi')
    plt.legend()
    plt.savefig('C:/Users/87588/Documents/RG Cosmo/Code/Figure/Phivsl%f.png'%c)
    plt.show()
    c=c+0.001

#%%
steps=[]
rArray=[]
tArray=[]
zArray=[]
phiArray=[]
rint=0.042
zint=0
tint=0
phiint=np.pi
flip=0
c=0.002
def evo(z,t,r,phi):
    global flip
    global c
    step=0.0001
    if(flip==0):
        if((z+1)**2-c**2/Y(r,t)**2<0):
            print('flip')
            flip=1
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
    return tempr,tempt,tempz,tempphi

for i in tqdm(range(0,240)):
    rArray.append(rint)
    tArray.append(tint)
    zArray.append(zint)
    phiArray.append(phiint)
    steps.append(i*0.01)
    rint, tint, zint, phiint=evo(zint,tint,rint,phiint)
    
string="c=%f"%c
plt.figure()
plt.plot(rArray,zArray,label=string) 
plt.xlabel('r')
plt.ylabel('z')
plt.legend()
plt.savefig('C:/Users/87588/Documents/RG Cosmo/Code/Figure/RvsZ%f.png'%c)
plt.show()

plt.figure()
plt.plot(steps,zArray,label=string) 
plt.xlabel('l')
plt.ylabel('z')
plt.legend()
plt.savefig('C:/Users/87588/Documents/RG Cosmo/Code/Figure/Zvsl%f.png'%c)
plt.show()

plt.figure()
plt.plot(steps,tArray,label=string) 
plt.xlabel('l')
plt.ylabel('t')
plt.legend()
plt.savefig('C:/Users/87588/Documents/RG Cosmo/Code/Figure/Tvsl%f.png'%c)
plt.show()

plt.figure()
plt.plot(steps,rArray,label=string) 
plt.xlabel('l')
plt.ylabel('r')
plt.legend()
plt.savefig('C:/Users/87588/Documents/RG Cosmo/Code/Figure/Rvsl%f.png'%c)
plt.show()

plt.figure()
plt.plot(steps,phiArray,label=string) 
plt.xlabel('l')
plt.ylabel('phi')
plt.legend()
plt.savefig('C:/Users/87588/Documents/RG Cosmo/Code/Figure/Phivsl%f.png'%c)
plt.show()
#%%
#F=[z,t,r,phi]
def f(l,F):
    c=0.1
    z=F[0]
    t=F[1]
    r=F[2]
    phi=F[3]
    f_z=-(dYdrdt(r,t)/dYdr(r,t))*((z+1)**2-c**2/Y(r,t)**2)-c**2*(dYdt(r,t)/Y(r,t)**3)
    f_t=z+1
    #f_r=(W(r)/dYdr(r,t))*np.sqrt((z+1)**2-c**2/Y(r,t)**2)
    f_r=fdr(r,t,z,c)
    f_phi=c**2/Y(r,t)**2
    return -1*np.array([f_z, f_t, f_r, f_phi])

zint=0
tint=0
rint=0.042
phiint=np.pi
y0=[zint,tint,rint,phiint]
l_span = np.linspace(0, 200, 201)
sol = solve_ivp(f, [0, 0.2], y0, method='Radau', l_eval=l_span)
plt.figure()
plt.plot(sol.y[2,:],sol.y[0,:]) 
plt.xlabel('r')
plt.ylabel('z')
plt.figure()
plt.plot(sol.t, sol.y[0,:],'-', label='z') 
plt.legend(loc='best')
plt.xlabel('l')
plt.figure()
plt.plot(sol.t, sol.y[2,:],'-', label='r')
plt.legend(loc='best')
plt.xlabel('l')
plt.figure()
plt.plot(sol.t, sol.y[1,:],'-', label='t')
plt.legend(loc='best')
plt.xlabel('l')
#%%
def f(l,F):
    c=0.1
    z=F[0]
    t=F[1]
    r=F[2]
    phi=F[3]
    f_z=-(dYdrdt(r,t)/dYdr(r,t))*((z+1)**2-c**2/Y(r,t)**2)-c**2*(dYdt(r,t)/Y(r,t)**3)
    f_t=z+1
    #f_r=(W(r)/dYdr(r,t))*np.sqrt((z+1)**2-c**2/Y(r,t)**2)
    f_r=fdr(r,t,z,c)
    f_phi=c**2/Y(r,t)**2
    return -1*np.array([f_z, f_t, f_r, f_phi])

zint=1.2
tint=-0.4
rint=0.042
phiint=np.pi
y0=[zint,tint,rint,phiint]
l_span = np.linspace(0, 200, 201)
sol = solve_ivp(f, [0, 0.1], y0, method='Radau', l_eval=l_span)
plt.figure()
plt.plot(sol.y[2,:],sol.y[0,:]) 
plt.xlabel('r')
plt.ylabel('z')
plt.figure()
plt.plot(sol.t, sol.y[0,:],'-', label='z') 
plt.legend(loc='best')
plt.xlabel('l')
plt.figure()
plt.plot(sol.t, sol.y[2,:],'-', label='r')
plt.legend(loc='best')
plt.xlabel('l')
plt.figure()
plt.plot(sol.t, sol.y[1,:],'-', label='t')
plt.legend(loc='best')
plt.xlabel('l')



#%%
t1=-0.2
r=np.linspace(0.0001,0.06,100)
y1=[]
dy1=[]
dy2=[]
dy3=[]
W0=[]
for i in range(len(r)):
    y1.append(Y(r[i],t1))
    dy1.append(dYdt(r[i],t1))
    dy2.append(dYdr(r[i],t1))
    dy3.append(dYdrdt(r[i],t1))
    W0.append(W(r[i]))
plt.figure()
plt.plot(r,y1,label='Y')
plt.legend(loc='best')
plt.figure()
plt.plot(r,dy1,label='dYdt')
plt.legend(loc='best')
plt.figure()
plt.plot(r,dy2,label='dYdr')
plt.legend(loc='best')
plt.figure()
plt.plot(r,dy3,label='dYdrdt')
plt.legend(loc='best')
plt.figure()
plt.plot(r,W0,label='W')
plt.legend(loc='best')
#%% plot diiferentials
t1=-0.2
z=0
c=0.01
r=np.linspace(-0.02,0.02,1000)
drdl=[]
root=[]
modDr=[]
for i in range(len(r)):
    drdl.append((W(r[i])/dYdr(r[i],t1))*np.sqrt((z+1)**2-c**2/Y(r[i],t1)**2))
    root.append((z+1)**2-c**2/Y(r[i],t1)**2)
    modDr.append(fdr(r[i],t1,z,c))
plt.figure()
plt.plot(r,drdl,label='drdl')
plt.legend(loc='best')
plt.figure()
plt.plot(r,root,label='root')
plt.legend(loc='best')
plt.figure()
plt.plot(r,modDr,label='modDr')
plt.legend(loc='best')
