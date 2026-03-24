# -*- coding: utf-8 -*-
"""
Created on Tue Mar 24 09:42:13 2026
@author: andrearaaschou
"""
import numpy as np
import  matplotlib.pyplot as plt

# Implement the deterministic and the stochastic (Gillespie algorithm) model for the bistable switch.
# Conduct simulations for both models with initial conditions specified on slides (21, 22, 23) 
# and with parameters specified in Gardner et al. 2000 Figure 5a (Reproduce the results shown in slides 21, 22, 23.)



# Initializing stuff i need but with random numbers - look up what they should be
a0 = 2 
a1 = 3
a2 = 4
a3 = 1
a4 = 4

ai = (a1, a2, a3, a4)

t0 = 4
x0 = 3

# Parameters from Gardner et al
alpha1 = 156.25 
alpha2 = 15.6

# more parameters from Gardner et al, not sure where to use these
beta = 2.5
gamma = 1
eta = 2.0015
K = 2.9618 * 10**(-5)

# ----------------- Start of program ---------------------
# Initialize t and X
t = t0
x = x0

# 1. Evaluate a_mu (a1, a2, a3, a4) and their sum a0 given the system is in state x at time t

# 2. Generate tao and mu using equation 8
# generate r1 and r2
r1, r2 = np.random.uniform(0, 1, 2)
print(f"r1: {r1}")
print(f"r2: {r2}")

tao = 1/a0 * np.log(1/r1)


#mu = np.min(ai[ai > r2 * a0]) # dont do this instead do
place = np.floor(r2 * a0) # round down to nearest integer
acumulative_conc = 0
for i in ai: # loop through ai
    acumulative_conc += ai[i] # add the current count
    if acumulative_conc >= place: # break the loop when the correct ai is found
        print(f'position is {i}')
        break
    
# take this element in the list or something similar

# 3. Update the system (t and x)

# 4. Record x, t. return to step 1 or else end simulation


# Results to reproduce:
    # one graph with the deterministic cell programming model 
    # one graph with the stochastic (Gillespie) model
    # time on x-axis, expression on y-axis
    # two different types of proteins expressed





















