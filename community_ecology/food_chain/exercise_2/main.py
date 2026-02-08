# -*- coding: utf-8 -*-
"""
Created on Sat Feb  7 10:55:29 2026
@author: andrearaaschou

Exercise 2 in food chain theory
"""
import numpy as np
import functions as f
from scipy.integrate import solve_ivp

# Model parameters
rR = 1      # Resource intrinsic growth rate
KR = 2      # Resource carrying kapacity
aC = 2      # Consumer attack rate
hC = 1.25   # Consumer handling time
eC = 1      # Conversion factor
muC = 0.4   # Consumer death rate
aP = 0.23   # Predator attack rate
hP = 2      # Predator handling time
eP = 1      # Predator conversion factor
muP = 0.1   # Predator death rate

# Initial conditions
R0 = 0.1
C0 = 0.1
P0 = 0.1
RCP0 = [R0, C0, P0]

t_interval = (0, 700) # time interval
t_eval = np.linspace(t_interval[0], t_interval[1], 1000) # density of points evaluated

# b) Run simulations of the system
out = solve_ivp(f.dRCPdt, 
                t_interval, 
                RCP0, 
                args=(rR, KR, aC, hC, eC, muC, aP, hP, eP, muP),
                t_eval = t_eval)

# Plot simulations
f.RC_simulation(out) 

# Print population sizes at the end of simulations
print("Equilibrium population sizes")
print(f" - Resource population size: {out.y.T[-1,0]:.2f}")
print(f" - Consumer population size: {out.y.T[-1,1]:.2f}")
print(f" - Predator population size: {out.y.T[-1,2]:.2f}")