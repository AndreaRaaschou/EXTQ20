# -*- coding: utf-8 -*-
"""
Created on Tue Mar 24 09:42:13 2026
@author: andrearaaschou

Exercise in gene regulatory networks:
Implement the deterministic and the stochastic (Gillespie algorithm) model for the bistable switch.
Conduct simulations for both models with initial conditions specified on slides (21, 22, 23) 
and with parameters specified in Gardner et al. 2000 Figure 5a
"""
import numpy as np
import  matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

def stochastic_simulation(t0, x0, alpha1, alpha2, beta, gamma, num_iterations):
    # empty arrays to store time and x
    t = np.zeros(num_iterations+1)
    x = np.zeros((num_iterations+1, 2))
    
    # ----------------- Start of program: Gillespie (stochastic) simulation ---------------------
    # 0. Initialize t and X
    t[0] = t0 
    x[0, :] = x0 
    
    for step in range(num_iterations):
        # 1. Evaluate a_mu (a1, a2, a3, a4) and their sum a0 given the system is in state x at time t
        a1 = alpha1 / (1 + x[step, 1]**beta)  # produce u
        a2 = x[step, 0]                       # degrade u
        a3 = alpha2 / (1 + x[step, 0]**gamma) # produce v
        a4 = x[step, 1]                       # degrade v
        a0 = a1 + a2 + a3 + a4 
        ai = [a1, a2, a3, a4]
        
        # 2. Generate tao and mu using equation 8
        r1, r2 = np.random.uniform(0, 1, 2) # generate r1 and r2
        tao = 1/a0 * np.log(1/r1) # generate tao
        
        # generate mu (get the position of the reaction that will occur)
        sum = 0
        position = -1
        for i in range(len(ai)):    # loop through ai
            sum += ai[i]            # add the current count
            if sum >= a0*r2:        # if the correct ai is found
                position = i        # save position
                break               # break the loop
        
        # 3. Update the system (t and x)
        t[step + 1] = t[step] + tao
        
        if position == 0:
            x[step + 1, 0] = x[step, 0] + 1
            x[step + 1, 1] = x[step, 1]
        elif position == 1:
            x[step + 1, 0] = x[step, 0] - 1
            x[step + 1, 1] = x[step, 1]
        elif position == 2:
            x[step + 1, 0] = x[step, 0]
            x[step + 1, 1] = x[step, 1] + 1 
        elif position == 3:
            x[step + 1, 0] = x[step, 0]
            x[step + 1, 1] = x[step, 1] - 1
            
    return (t, x)

def gillespie_model(t, x, alpha1, alpha2, beta, gamma):
    u, v = x
    dudt = alpha1 / (1 + v**beta) - u 
    dvdt = alpha2 / (1 + u**gamma) - v
    return (dudt, dvdt)

def make_plots(t, x, sol):
    fig, (ax1, ax2) = plt.subplots(1,2)
    
    ax1.plot(sol.t, sol.y[0], label = 'u')
    ax1.plot(sol.t, sol.y[1], label = 'v')
    ax1.set_title('Deterministic Gillespie model', fontsize=10)
    ax1.set_xlabel('Time', fontsize=10)
    ax1.set_ylabel('number of molecules', fontsize=10)
    ax1.legend(loc='upper right')
    
    ax2.plot(t, x[:, 0], label = 'u')
    ax2.plot(t, x[:, 1], label = 'v')
    ax2.set_title('Stochastic Gillespie simulation', fontsize=10)
    ax2.set_xlabel('Time', fontsize=10)
    ax2.set_ylabel('number of molecules', fontsize=10)
    ax2.legend(loc='upper left')
    
    fig.tight_layout(pad=3.0)
    return (fig, ax1, ax2)

def run_simulations(t0, x0):
    # Parameters from Gardner et al
    alpha1 = 156.25    # effective rate of synthesis of repressor 1
    alpha2 = 15.6      # effective rate of synthesis of repressor 2
    beta = 2.5         # cooperativity of repression of promoter 2
    gamma = 1          # cooperativity of repression of promoter 1
    
    # deterministic model
    t_span = (0, 20)
    t_eval = np.linspace(*t_span, 100)
    sol = solve_ivp(gillespie_model, t_span, x0, t_eval=t_eval, args = (alpha1, alpha2, beta, gamma))
    
    # stochastic model
    t, x = stochastic_simulation(t0, x0, alpha1, alpha2, beta, gamma, num_iterations = 200)
    
    # make plots
    fig, ax1, ax2 = make_plots(t, x, sol)
    fig.show()

# run the simulation with two different initial conditions
run_simulations(0, np.array((6, 1)))
run_simulations(0, np.array((1, 6)))


















