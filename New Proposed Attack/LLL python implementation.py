# -*- coding: utf-8 -*-
"""
Created on Sat Oct  2 21:25:39 2021

@author: Charissa
"""

#using lll attack

import numpy as np
import math

l1 = np.array([1,0,1,0], dtype = ('float64'))
l2 = np.array([1,1,1,1], dtype = ('float64'))
l3 = np.array([0,1,2,1],dtype = ('float64'))
latticelist = [l1,l2,l3]


def norm(vector):
    num = 0
    for i in range(len(vector)):
        num += (vector[i])**2
    normvec = num**(1/2)
    return normvec



def gramschmidts(latticelist):
    '''
    

    Parameters
    ----------
    latticelist : list
        a list of basis of lattice point of NTRU 

    Returns
    -------
    ortholist : list
        outputs a list of orthonormal basis

    '''
    newlist = []
    vector1 = latticelist[0]
    newlist.append(vector1)
    for i in range(1,len(latticelist)):
        w= latticelist[i]
        nextvector =0
        nextvector += latticelist[i]
        for j in range(len(newlist)):
            nextvector -= (np.dot(w,newlist[j])/((norm(newlist[j]))**2)) * newlist[j]
        newlist.append(nextvector)
    ortholist = []
    for k in range(len(newlist)):
        ortholist.append(newlist[k]/norm(newlist[k]))
    return ortholist


def lovasz(ortholist):
    boolean = False
    for i in range(len(ortholist)):
        for j in range(len(ortholist)):
            if i>j:
                if np.dot(ortholist[i],ortholist[j])/((norm(ortholist[j]))**2) <= 0.5 :
                    pass
                else:
                    temp = ortholist[i]
                    ortholist[i] = ortholist[j]
                    ortholist[j] = temp                
    for k in range(1,len(ortholist)):
        if (norm(ortholist[k-1]))**2> (norm(ortholist[k]))**2 + ((np.dot(ortholist[k],ortholist[k-1])/((norm(ortholist[k-1]))**2))**2)*norm((ortholist[k-1]))**2:
            temp = ortholist[k-1]
            ortholist[k-1] = ortholist[k]
            ortholist[k] = temp 
            boolean = True
    return ortholist, boolean

def lllalgo(l):
   if lovasz(gramschmidts(l))[1] == True:
       lllalgo(lovasz(gramschmidts(l))[0])
   else:
        return lovasz(gramschmidts(l))[0]
        
print(lllalgo(latticelist))   
     

    