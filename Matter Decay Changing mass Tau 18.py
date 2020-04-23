#!/usr/bin/env python
# coding: utf-8

# In[1]:


import numpy as np
import math
from scipy.integrate import quad
from array import array
import matplotlib.pyplot as plt
from sympy.solvers import solve
from sympy import *
from scipy.integrate import odeint
from IPython.core.debugger import set_trace


# In[30]:


#parameters
#H_0 = 2.26853E-18 #seconds-1
H_0 = 0.07152 # gyr-1
c = 9.7156E-15 #Mpc/s
p_crit = 9.47E-27 #kg/m^3
p_m0 = 0.2722*p_crit
tau = 18 #gyr (100 billion yrs)
omega_m = 0.1684 #Best overal fit 0.0863
omega_l = 0.7228
omega_r = 8.4E-5


# In[31]:


def model(z,t):
    x = z[0]
    y = z[1]
    if z[0]== 0 or z[1] == 0:
        dxdt = 3/2*H_0*math.sqrt(omega_m*exp((13.87-t)/tau)+omega_l*x**2)
        dydt = 0
    else:
        dxdt = 3/2*H_0*math.sqrt(omega_m*exp((13.87-t)/tau)+omega_l*x**2+y/(x**(2/3)))
        dydt = omega_m*(x**(2/3))*exp((13.87-t)/tau)/tau
    dzdt = [dxdt,dydt]
    return dzdt


# In[32]:


# initial condition
z0 = [0,0]

# number of time points
n = 10000

# time points
t = np.linspace(0,14,n)


# store solution
#x = np.empty_like(t)
#y = np.empty_like(t)
x = np.zeros(n)
y = np.zeros(n)
# record initial conditions
x[0] = z0[0]
y[0] = z0[1]


# In[33]:


# solve ODE
for i in range(1,n):
    # span for next time step
    tspan = [t[i-1],t[i]]
    
    # solve for next step
    z = odeint(model,z0,tspan)
    #set_trace()
    if z[1][0]<0:
        break
    # store solution for plotting
    x[i] = z[1][0]
    y[i] = z[1][1]
    if z[1][0]==0:
        break
    # next initial condition
    z0 = z[1]


# In[34]:


a = x**(2/3)


# In[35]:


get_ipython().run_line_magic('store', '-r t_stand')
get_ipython().run_line_magic('store', '-r a_stand')


# In[36]:


fig = plt.figure()
ax1 = fig.add_subplot(111)
ax2 = fig.add_subplot(111)
#ax2 = ax1.twiny()
#ax1.set_xlim(-15,1)
#ax1.set_ylim(0,2)
ax1.plot(t,a, label= 'Tau 10')
ax2.plot(t_stand,a_stand,linestyle = '--', label = 'Standard Model')
#ax1.set_xticks([-1,-0.5,0,0.5,1])
#ax1.set_xticklabels([-1,-0.5,0,0.5,1])
ax1.set_xlabel(r"Time (Gy)")
ax1.set_ylabel(r"Scale Factor, a")
#ax1Xs = ax1.get_xticks()
#ax2Xs = []
#for X in ax1Xs:
    #ax2Xs.append(X * 13.75)
#ax2.set_xticks(ax1Xs)
#ax2.set_xbound(ax1.get_xbound())
#ax2.set_xticklabels(ax2Xs)
#ax2.set_xlabel(r"Look-Back Time (Gy)")
plt.legend()
plt.show()


# In[37]:


k = abs(a_stand[1:] - a[1:])


# In[38]:


test = k/a_stand[1:]


# In[39]:


if np.amax(test)<= 0.1:
    print('pass')
else:
    print('Fail')


# In[40]:


np.amax(test)


# In[41]:


np.where(test == np.amax(test))


# In[42]:


a_18 = a
get_ipython().run_line_magic('store', 'a_18')


# In[ ]:




