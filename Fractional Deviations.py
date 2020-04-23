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


# In[2]:


#parameters
#H_0 = 2.26853E-18 #seconds-1
H_0 = 0.07152 # gyr-1
c = 9.7156E-15 #Mpc/s
p_crit = 9.47E-27 #kg/m^3
p_m0 = 0.2722*p_crit 
omega_l = 0.7228
omega_r = 8.4E-5


# In[16]:


get_ipython().run_line_magic('store', '-r a_stand')
get_ipython().run_line_magic('store', '-r t_stand')
get_ipython().run_line_magic('store', '-r a_10')
get_ipython().run_line_magic('store', '-r a_20')
get_ipython().run_line_magic('store', '-r a_100')
get_ipython().run_line_magic('store', '-r a_1000')


# In[17]:


a10 = a_10/a_stand
a10[0]= 1
a20 = a_20/a_stand
a20[0] = 1
a100 = a_100/a_stand
a100[0]= 1
a1000 = a_1000/a_stand
a1000[0]=1


# In[18]:


fig = plt.figure()
ax1 = fig.add_subplot(111)
ax2 = fig.add_subplot(111)
ax3 = fig.add_subplot(111)
ax4 = fig.add_subplot(111)
#ax2 = ax1.twiny()
#ax1.set_xlim(-15,1)
#ax1.set_ylim(0,2)
ax1.plot(t_stand,a10, label = 'Tau 10')
ax2.plot(t_stand,a20, label= 'Tau 20')
ax3.plot(t_stand, a100, label = 'Tau 100')
ax4.plot(t_stand, a1000, label = 'Tau 1000')
ax1.set_xlabel(r"Time (Gy)")
ax1.set_ylabel(r"Scale Factor Ratio (Decay/Standard)")
plt.legend()
plt.show()


# In[ ]:




