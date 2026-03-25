# -*- coding: utf-8 -*-
"""
Created on Wed Mar 25 16:05:16 2026
@author: andrearaaschou
"""
import numpy as np
import  matplotlib.pyplot as plt

def stochastic_simulation(t0, x0, alpha1, alpha2, beta, gamma, num_iterations):
    # empty arrays to store time and x
    t = np.zeros(num_iterations+1)
    x = np.zeros((num_iterations+1, 2))
    
    # ----------------- Start of program: Gillespie (stochastic) simulation ---------------------
    # 0. Initialize t and X
    t[0] = t0 
    x[0, :] = x0 
    
    for step in range(num_iterations):
        #print("\nNEW LOOP")
        
        # 1. Evaluate a_mu (a1, a2, a3, a4) and their sum a0 given the system is in state x at time t
        a1 = alpha1 / (1 + x[step, 1]**beta)  # produce u
        a2 = x[step, 0]                       # degrade u
        a3 = alpha2 / (1 + x[step, 0]**gamma) # produce v
        a4 = x[step, 1]                       # degrade v
        a0 = a1 + a2 + a3 + a4 
        ai = [a1, a2, a3, a4]
        
        #print(f"a1: {a1}, a2: {a2}, a3: {a3}, a4: {a4}")
        #print(f"a0: {a0}")
        
        
        # 2. Generate tao and mu using equation 8
        r1, r2 = np.random.uniform(0, 1, 2) # generate r1 and r2
        #print(f"r1: {r1}")
        #print(f"r2: {r2}")
        tao = 1/a0 * np.log(1/r1) # generate tao
        
        # generate mu (get the position of the reaction that will occur)
        sum = 0
        position = -1
        for i in range(len(ai)):    # loop through ai
            sum += ai[i]            # add the current count
            #print(f"sum: {sum}")
            if sum >= a0*r2:        # if the correct ai is found
                position = i        # save position
                #print(f'position is {i}')
                break               # break the loop
        
        
        
        # 3. Update the system (t and x)
        t[step + 1] = t[step] + tao
        #print(f"t: {t}")
        
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
        #print(f"x: {x}")    
            
    return (t, x)
