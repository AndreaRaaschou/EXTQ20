# -*- coding: utf-8 -*-
"""
Created on Wed Feb 18 10:48:58 2026
@author: andrearaaschou
"""
import numpy as np
import  matplotlib.pyplot as plt

# Task 1
# Make plots of the given data in strokework

t = np.arange(10)
 

plt.plot(t, np.sin(t))
   
plt.title('matplotlib.pyplot.ginput()\
 function Example', fontweight ="bold")
 
print("After 3 clicks :")
x = plt.ginput(3)
print(x)
 
plt.show()
