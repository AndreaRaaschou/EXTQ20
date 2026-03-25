# -*- coding: utf-8 -*-
"""
Created on Tue Mar 24 09:42:13 2026
@author: andrearaaschou

Implement the deterministic and the stochastic (Gillespie algorithm) model for the bistable switch.
Conduct simulations for both models with initial conditions specified on slides (21, 22, 23) 
and with parameters specified in Gardner et al. 2000 Figure 5a (Reproduce the results shown in slides 21, 22, 23.)
"""
import numpy as np
import  matplotlib.pyplot as plt
import gillespie_simulation as gs

# Parameters from Gardner et al
alpha1 = 156.25    # effective rate of synthesis of repressor 1
alpha2 = 15.6      # effective rate of synthesis of repressor 2
beta = 2.5         # cooperativity of repression of promoter 2
gamma = 1          # cooperativity of repression of promoter 1
  
# more parameters from Gardner et al, not sure where to use these - will probably remove?
#eta = 2.0015
#K = 2.9618 * 10**(-5)

# Initial conditions to test and put in report:
    # U(0) = 6, V(0) = 1
    # U(0) = 1, V(0) = 6
t0 = 0                  # starting time of the simulation
x0 = np.array((6, 1))   # initial number of molecules (u, v)

def gillespie_model(t, x, alpha1, alpha2, beta, gamma):
    u, v = x
    dudt = alpha1 / (1 + v**beta) - u 
    dvdt = alpha2 / (1 + u**gamma) - v
    return (dudt, dvdt)



t, x = gs.stochastic_simulation(t0, x0, alpha1, alpha2, beta, gamma, num_iterations = 200)

fig, (ax1, ax2) = plt.subplots(1,2)
ax2.plot(t, x[:, 0], label = 'u')
ax2.plot(t, x[:, 1], label = 'v')
ax2.set_title('plot title', fontsize=12)
ax2.set_xlabel('Time', fontsize=10)
ax2.set_ylabel('number of molecules', fontsize=10)
ax2.legend()























