# -*- coding: utf-8 -*-
"""
Created on Sat Feb  7 11:16:59 2026
@author: andrearaaschou
"""
import  matplotlib.pyplot as plt
    
#' b) Calculates the growth rates of the three populations
def dRCPdt(t, RCPvector, rR, KR, aC, hC, eC, muC, aP, hP, eP, muP):
    R = RCPvector[0] # resource density
    C = RCPvector[1] # consumer density
    P = RCPvector[2] # predator density
    
    dRdt = rR*R * (1 - R/KR) - aC*R*C/(1 + aC*hC*R)
    dCdt = eC*aC*R*C/(1 + aC*hC*R) - muC*C - aP*C*P/(1 + aP*hP*C)
    dPdt = eP*aP*C*P/(1 + aP*hP*C) - muP*P
    
    return [dRdt, dCdt, dPdt]

def RC_simulation(out):
    plt.plot(out.t, out.y.T[:,0], label="Resource")
    plt.plot(out.t, out.y.T[:,1], label = "Consumer")
    plt.plot(out.t, out.y.T[:,2], label = "Predator")
    plt.xlabel("time")
    plt.ylabel("population size")
    plt.legend(loc="upper right")
    plt.show()
