# -*- coding: utf-8 -*-
"""
Created on Sat Oct  2 21:25:39 2021

@author: Charissa
"""

#using lll attack

import numpy as np
from numpy.linalg import norm 
import math


l1 = np.array([1,0,0,0,0,0,0,0,0,0,0,50,5,32,36,31,53,28,46,25,49,11])
l2 = np.array([0,1,0,0,0,0,0,0,0,0,0,11,50,5,32,36,31,53,28,46,25,49])
l3 = np.array([0,0,1,0,0,0,0,0,0,0,0,49,11,50,5,32,36,31,53,28,46,25])
l4 = np.array([0,0,0,1,0,0,0,0,0,0,0,25,49,11,50,5,32,36,31,53,28,46])
l5 = np.array([0,0,0,0,1,0,0,0,0,0,0,46,25,49,11,50,5,32,36,31,53,28])
l6 = np.array([0,0,0,0,0,1,0,0,0,0,0,28,46,25,49,11,50,5,32,36,31,53])
l7 = np.array([0,0,0,0,0,0,1,0,0,0,0,53,28,46,25,49,11,50,5,32,36,31])
l8 = np.array([0,0,0,0,0,0,0,1,0,0,0,31,53,28,46,25,49,11,50,5,32,36])
l9 = np.array([0,0,0,0,0,0,0,0,1,0,0,36,31,53,28,46,25,49,11,50,5,32])
l10 = np.array([0,0,0,0,0,0,0,0,0,1,0,32,36,31,53,28,46,25,49,11,50,5])
l11 = np.array([0,0,0,0,0,0,0,0,0,0,1,5,32,36,31,53,28,46,25,49,11,50])
l12 = np.array([0,0,0,0,0,0,0,0,0,0,0,61,0,0,0,0,0,0,0,0,0,0])
l13 = np.array([0,0,0,0,0,0,0,0,0,0,0,0,61,0,0,0,0,0,0,0,0,0])
l14 = np.array([0,0,0,0,0,0,0,0,0,0,0,0,0,61,0,0,0,0,0,0,0,0])
l15 = np.array([0,0,0,0,0,0,0,0,0,0,0,0,0,0,61,0,0,0,0,0,0,0])
l16 = np.array([0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,61,0,0,0,0,0,0])
l17 = np.array([0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,61,0,0,0,0,0])
l18 = np.array([0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,61,0,0,0,0])
l19 = np.array([0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,61,0,0,0])
l20 = np.array([0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,61,0,0])
l21 = np.array([0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,61,0])
l22 = np.array([0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,61])


latticelist= [l1,l2,l3,l4,l5,l6,l7,l8,l9,l10,l11,l12,l13,l14,l15,l16,l17,l18,l19,l20,l21,l22]

'''
l1 = np.array([1,1,1])
l2 = np.array([-1,0,2])
l3 = np.array([3,5,6])
latticelist= [l1,l2,l3]

'''
#finding norm
# use norm(vector,2) for norm instead from numpy.linalg



def gramschmidts(latticelist):
    '''
    

    Parameters
    ----------
    latticelist : list
        a list of basis of lattice point of NTRU 

    Returns
    -------
    ortholist : list
        outputs a list of orthogonal basis

    '''
    newlist = []
    newlist.append(latticelist[0])
    for i in range(1,len(latticelist)):
        nextvector = np.zeros(len(latticelist[i]))
        nextvector = latticelist[i]
        for j in range(len(newlist)):
            nextvector = nextvector - round(np.dot(latticelist[i],newlist[j])/(norm(newlist[j],2))**2) * newlist[j]
        newlist.append(nextvector)
    return newlist


def lllalgo(l):
    basis = gramschmidts(l)
    for j in range(len(l)+1):
        for k in range(j+1,len(l)):
            coef = np.dot(l[k],basis[j])/(norm(basis[j],2))**2
            if abs(round(coef)) > 0.5:
                l[k] = l[k] - round(coef)*l[j]
                basis = gramschmidts(l)
    for i in range(len(l)-1):
        if 0.75*norm(basis[i],2)**2> norm(basis[i+1]+(np.dot(l[i],basis[i+1])/(norm(basis[i+1],2))**2)*basis[i],2)**2:
            temp = l[i]
            l[i]= l[i+1]
            l[i+1]= temp
            lllalgo(l)
    else:
        return l

print(lllalgo(latticelist))

print()

largestnorm =[]
for i in range(len(lllalgo(latticelist))):
    largestnorm.append(norm(lllalgo(latticelist)[i],2))
print(max(largestnorm))
print(largestnorm)

