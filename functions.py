#!/usr/bin/env python
# coding: utf-8

# In[1]:


import numpy as np
import math
from scipy.integrate import quad
from array import array
from sympy.solvers import solve
from sympy import *
from scipy.integrate import odeint
from IPython.core.debugger import set_trace
from scipy.optimize import fsolve


# In[7]:


def omega_find(tau,m_min,m_max):
    m = 20
    array_m = np.linspace(m_min,m_max,m)
    results = np.zeros((m,2))
    
    for q in range(1,m+1):
        omega_m = array_m[q-1]


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
    
        #Solve ODE
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
    
        a = x**(2/3)
    
        k = abs(a_stand[1:] - a[1:])
    
        test = k/a_stand[1:]
        
        if np.amax(test)<= 0.1:
            results[q-1][0] = np.amax(test)
            results[q-1][1] = omega_m
        else:
            results[q-1][0] = 1
            results[q-1][1] = omega_m 
        
        min_err = np.amin(results[:-1],axis =0)[0]
    
        index = np.where(results == min_err)[0]
    
        omega_m = results[index][0][1]
    return omega_m


# In[8]:


#parameters
#H_0 = 2.26853E-18 #seconds-1
H_0 = 0.07152 # gyr-1
c = 9.7156E-15 #Mpc/s
p_crit = 9.47E-27 #kg/m^3
p_m0 = 0.2722*p_crit
tau = 10 #gyr (100 billion yrs)
omega_l = 0.7228
omega_r = 8.4E-5


# In[9]:


get_ipython().run_line_magic('store', '-r t_stand')
get_ipython().run_line_magic('store', '-r a_stand')


# In[ ]:


#Find Best Omega_m,0 Value in interval
omega_find(tau,0.1,0.2)


# ## FInding Omega_m using shooting method

# In[15]:


#parameters
#H_0 = 2.26853E-18 #seconds-1
H_0 = 0.07152 # gyr-1
c = 9.7156E-15 #Mpc/s
p_crit = 9.47E-27 #kg/m^3
p_m0 = 0.2722*p_crit
tau = 18 #gyr (100 billion yrs)
omega_l = 0.7228
omega_r = 8.4E-5


# In[16]:



def model(z,t):
            x,y,omega_m = z
            if z[0]== 0 or z[1] == 0:
                dxdt = 3/2*H_0*math.sqrt(omega_m*exp((13.87-t)/tau)+omega_l*x**2)
                dydt = 0
            else:
                dxdt = 3/2*H_0*math.sqrt(omega_m*exp((13.87-t)/tau)+omega_l*x**2+y/(x**(2/3)))
                dydt = omega_m*(x**(2/3))*exp((13.87-t)/tau)/tau
            dzdt = [dxdt,dydt,0]
            return dzdt
# initial condition
z0 = [0,0]
# number of time points
n = 1000

# time points
t = np.linspace(0,13.87,n)


# In[17]:


def objective(omega_m):
    z0 = [0,0,omega_m]
    tspan = np.linspace(0, 13.87,1000)
    z = odeint(model, z0, tspan)
    x = z[:,0]
    return x[-1] -1


# In[18]:


omega_m = fsolve(objective, 0.1)


# In[19]:


print(omega_m)


# In[ ]:




