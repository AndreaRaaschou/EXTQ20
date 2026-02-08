# -*- coding: utf-8 -*-
"""
Created on Sat Feb  7 08:42:27 2026
@author: andrearaaschou
"""
import numpy as np
import functions as f
from scipy.integrate import solve_ivp

# Model parameters
rR = 1      # Resource intrinsic growth rate
KR = 1      # Resource carrying kapacity
aC = 2      # Consumer attack rate
hC = 1.25   # Consumer handling time
eC = 1      # Conversion factor
muC = 0.4   # Consumer death rate

# Initial conditions
R0 = 0.1
C0 = 0.1
RC0 = [R0, C0]

t_interval = (0, 100) # time interval
t_eval = np.linspace(t_interval[0], t_interval[1], 1000) # density of points evaluated

# b) Run simulations of the system
out = solve_ivp(f.dRCdt, 
                t_interval, 
                RC0, 
                args=(rR, KR, aC, hC, eC, muC),
                t_eval = t_eval)

# Plot isoclines
f.RC_isoclines(rR, KR, aC, hC, eC, muC, out)

# Plot simulations
f.RC_simulation(out) 

# Print population sizes at the end of simulations
print("Equilibrium population sizes")
print(f" - Resource population size: {out.y.T[-1,0]:.2f}")
print(f" - Consumer population size: {out.y.T[-1,1]:.2f}")

