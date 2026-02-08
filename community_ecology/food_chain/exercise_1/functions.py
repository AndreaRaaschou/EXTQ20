# -*- coding: utf-8 -*-
"""
Created on Sat Feb  7 08:39:14 2026
@author: andrearaaschou
"""
import numpy as np
import  matplotlib.pyplot as plt

#' a) Plots the isoclines of a short food chain of one resource
#' Isoclines have been solved analytically
def RC_isoclines(rR, KR, aC, hC, eC, muC, out):
    # Resource isocline (R, C(R))
    RR = np.linspace(0, KR, 100)  # Resource-axis (0, carrying kapacity)
    CR = (rR/aC) * (1 - RR/KR) * (1 + aC*hC*RR)
    
    # Consumer isocline (verticle line)
    RC = muC / (eC*aC - muC*aC*hC)
    
    # plot isoclines
    plt.plot(RR, CR, label="Resource isocline")
    plt.axvline(RC, color='red', label="Consumer isocline")
    plt.xlabel("R")
    plt.ylabel("C")
    
    # Trajectory
    plt.plot(out.y.T[:,0], out.y.T[:,1], 'k', label="Trajectory")
    plt.plot(out.y.T[0,0], out.y.T[0,1], '*k', label="Starting point")
    
    plt.legend(loc="upper right")
    plt.show()
    
#' b) Calculates the growth rates of the two populations
def dRCdt(t, RCvector, rR, KR, aC, hC, eC, muC):
    R = RCvector[0] # resource density
    C = RCvector[1] # consumer density
    
    dRdt = rR*R * (1 - R/KR) - aC*R*C/(1 + aC*hC*R)
    dCdt = eC*aC*R*C/(1 + aC*hC*R) - muC*C
    
    return [dRdt, dCdt]

def RC_simulation(out):
    plt.plot(out.t, out.y.T[:,0], label="Resource")
    plt.plot(out.t, out.y.T[:,1], label = "Consumer")
    plt.xlabel("time")
    plt.ylabel("population size")
    plt.legend(loc="upper right")
    plt.show()
