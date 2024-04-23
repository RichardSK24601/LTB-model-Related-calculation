import math
import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt

# y = np.array([r, u, phi, time])
def f(t, y):
    r = y[0] 
    f_r = y[1] # this is the dr / dT auxiliary equation
    f_u = - 7.5 / (r**2) + 2 / (r**3) - 15 / (r**4)
    f_phi = 1 / (r**2)
    f_time = k1 * r / (r - 5) # this is the equation of the time coordinate
    return np.array([f_r, f_u, f_phi, f_time])

# from the initial value for r = r0 and given energy k,  
# calculate the initial rate of change dr / dT = u0
def ivp(r0, k, sign):
    u0 = math.sqrt( k - ( 1 - 5 / (r0**2) ) * ( 3 + 2 / (r0**2) ) )
    return sign * u0

k = 3.0
k1 = 2.0
r0 = 20.0
sign = 1 # or -1

u0 = ivp(r0, k, sign)
# y = np.array([r, u, phi, time])
y0 = [r0, u0, math.pi, 0]
t_span = np.linspace(0, 1000, num=1001)   

sol = solve_ivp(f, [0, 1000], y0, method='Radau', t_eval=t_span)

plt.plot(sol.t, sol.y[0,:],'-', label='r(t)') 
plt.plot(sol.t, sol.y[2,:],'-', label='phi(t)')
plt.legend(loc='best')
plt.xlabel('T')

#%%