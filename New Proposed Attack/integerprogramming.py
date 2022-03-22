# -*- coding: utf-8 -*-
"""
Created on Mon Mar 21 21:43:08 2022

@author: Charissa
"""

"Linear programming"

from scipy import optimize
import numpy as np
import matplotlib.pyplot as plt



c = np.array([0, 0, 0, 0, 0, -1, -1, -1, -1, -1, -1, -1])
A_ub = np.array([
        [0, 0, -20, -14, 0, 1, 0, 0, 0, 0, 0, 0],
        [-1, -6, 0, 0, -15, 0, 1, 0, 0, 0, 0, 0],
        [-2, -12, 0, 0, -13, 0, 0, 1, 0, 0, 0, 0],
        [-3, -1, 0, 0, -2, 0, 0, 0, 1, 0, 0, 0],
        [-3, -5, 0, 0, -30, 0, 0, 0, 0, 1, 0, 0],
        [0, 0, -8, -3, -30, 0, 0, 0, 0, 0, 1, 0],
        [-10, 0, -30, -50, 0, 0, 0, 0, 0, 0, 0, 1]])

A_eq = np.array([[1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0]])
b_ub = np.array([-11, -14, -16, -40, -60, -12, -12])

bounds = [(0, None), (0, None), (0, None), (0, None), (0, None), (0, 29), (0, 36), (0, 29), (0, 50), (0, 30), (0, 128), (0, 78)]



def lp(c,A_ub,A_eq,b_ub,bounds):
    for i in range(16, 20):                    
        b_eq = np.array([i])
        print("For value of " + str(i)+', the ')
        a = optimize.linprog(c, A_ub = A_ub, b_ub = b_ub, A_eq = A_eq, b_eq = b_eq, bounds = bounds)

        print(a.fun)
        print()

    # print((optimize.linprog(c, A_ub = A_ub, b_ub = b_ub, A_eq = A_eq, b_eq = np.array([11]), bounds = bounds)))    
    # print((optimize.linprog(c, A_ub = A_ub, b_ub = b_ub, A_eq = A_eq, b_eq = np.array([11.0]), bounds = bounds)))    

    x = np.linspace(16, 19, 1000)   
    y = []
    for item in range(len(x)):
        y.append((optimize.linprog(c, A_ub = A_ub, b_ub = b_ub, A_eq = A_eq, b_eq = np.array([x[item]]), bounds = bounds)).fun)
    # y = (optimize.linprog(c, A_ub = A_ub, b_ub = b_ub, A_eq = A_eq, b_eq = np.array([x]), bounds = bounds)).fun
    
    f = plt.figure(figsize=(16, 8))
    plt.plot(x,y)
    plt.show()



lp(c,A_ub,A_eq,b_ub,bounds)