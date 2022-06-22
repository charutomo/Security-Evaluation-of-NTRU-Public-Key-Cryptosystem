# -*- coding: utf-8 -*-
"""
Created on Mon Mar 21 21:43:08 2022

@author: Charissa
"""


#backtracking 
import numpy as np
import math 
import copy
from numpy.linalg import norm 
from conversiontxtfile import converttolist
from conversiontxtfile import converttoarray
from checkdata import SVFinder


latlist = converttolist('testdata.txt')
def gramschmidts(latticelist):
    '''
    gram schimdts orthogonalisation

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
initial = SVFinder(gramschmidts(latlist))[0]


def checknorm(l):
    val = np.zeros(22)
    for i in range(len(l)):
        val+=l[i]*latlist[i]
    return norm(val,2)

for i in range(len(initial)):
    temp = 0
    left = math.floor(initial[i])
    leftlist = copy.deepcopy(initial)
    leftlist[i] = left
    right = round(left+1)
    rightlist = copy.deepcopy(initial)
    rightlist[i] = right
    if checknorm(leftlist)<checknorm(rightlist):
        initial = leftlist
    elif checknorm(leftlist)>checknorm(rightlist):
        initial = rightlist
    else:
        initial = leftlist 
        
print(initial, checknorm(initial))    







