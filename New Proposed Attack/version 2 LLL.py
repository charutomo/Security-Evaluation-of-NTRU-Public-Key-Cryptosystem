# -*- coding: utf-8 -*-
"""
Created on Wed Mar  9 13:48:42 2022

@author: Charissa
"""

import numpy as np
from numpy.linalg import norm 
import math
from conversiontxtfile import convert 
from checkdata import SV

latticelist = convert("testdata.txt")


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
    k = 1
    while k<= len(l)-1:
        for j in range(k-1,-1,-1):
            coef = np.dot(l[k],basis[j])/(norm(basis[j],2))**2
            if abs(coef) > 0.5:
                l[k] = l[k] -  round(coef)* l[j]
                basis = gramschmidts(l)
        if norm(basis[k],2)**2 >= (0.75-(np.dot(l[k],basis[k-1])/(norm(l[k-1],2)**2))**2)*norm(basis[k-1],2)**2:
            k = k+1
        else:
            temp = l[k]
            l[k]= l[k-1]
            l[k-1]= temp
            k = max(k-1,1)
    return l


print(lllalgo(latticelist))


print()

#check for largest norm which determines the SVP
largestnorm =[]
for i in range(len(lllalgo(latticelist))):
    largestnorm.append(norm((lllalgo(latticelist))[i],2))
print(max(largestnorm))
print(largestnorm)
