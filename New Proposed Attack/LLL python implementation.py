# -*- coding: utf-8 -*-
"""
Created on Sat Oct  2 21:25:39 2021

@author: Charissa
"""

#using lll attack

import numpy as np
from numpy.linalg import norm 
import math

l1 = np.array([1,0,0,0,0,0,0,0,0,0,0,50,5,32,36,31,53,28,46,25,49,11], dtype = ('float64'))
l2 = np.array([0,1,0,0,0,0,0,0,0,0,0,11,50,5,32,36,31,53,28,46,25,49], dtype = ('float64'))
l3 = np.array([0,0,1,0,0,0,0,0,0,0,0,49,11,50,5,32,36,31,53,28,46,25], dtype = ('float64'))
l4 = np.array([0,0,0,1,0,0,0,0,0,0,0,25,49,11,50,5,32,36,31,53,28,46], dtype = ('float64'))
l5 = np.array([0,0,0,0,1,0,0,0,0,0,0,46,25,49,11,50,5,32,36,31,53,28], dtype = ('float64'))
l6 = np.array([0,0,0,0,0,1,0,0,0,0,0,28,46,25,49,11,50,5,32,36,31,53], dtype = ('float64'))
l7 = np.array([0,0,0,0,0,0,1,0,0,0,0,53,28,46,25,49,11,50,5,32,36,31], dtype = ('float64'))
l8 = np.array([0,0,0,0,0,0,0,1,0,0,0,31,53,28,46,25,49,11,50,5,32,36], dtype = ('float64'))
l9 = np.array([0,0,0,0,0,0,0,0,1,0,0,36,31,53,28,46,25,49,11,50,5,32], dtype = ('float64'))
l10 = np.array([0,0,0,0,0,0,0,0,0,1,0,32,36,31,53,28,46,25,49,11,50,5], dtype = ('float64'))
l11 = np.array([0,0,0,0,0,0,0,0,0,0,1,5,32,36,31,53,28,46,25,49,11,50], dtype = ('float64'))
l12 = np.array([0,0,0,0,0,0,0,0,0,0,0,61,0,0,0,0,0,0,0,0,0,0], dtype = ('float64'))
l13 = np.array([0,0,0,0,0,0,0,0,0,0,0,0,61,0,0,0,0,0,0,0,0,0], dtype = ('float64'))
l14 = np.array([0,0,0,0,0,0,0,0,0,0,0,0,0,61,0,0,0,0,0,0,0,0], dtype = ('float64'))
l15 = np.array([0,0,0,0,0,0,0,0,0,0,0,0,0,0,61,0,0,0,0,0,0,0], dtype = ('float64'))
l16 = np.array([0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,61,0,0,0,0,0,0], dtype = ('float64'))
l17 = np.array([0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,61,0,0,0,0,0], dtype = ('float64'))
l18 = np.array([0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,61,0,0,0,0], dtype = ('float64'))
l19 = np.array([0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,61,0,0,0], dtype = ('float64'))
l20 = np.array([0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,61,0,0], dtype = ('float64'))
l21 = np.array([0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,61,0], dtype = ('float64'))
l22 = np.array([0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,61], dtype = ('float64'))

latticelist= [l1,l2,l3,l4,l5,l6,l7,l8,l9,l10,l11,l12,l13,l14,l15,l16,l17,l18,l19,l20,l21,l22]



#finding norm
# use norm(vector,2) for norm instead from numpy.linalg
# for norm^2 use norm(vector,1)


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
        w= latticelist[i]
        nextvector =0
        nextvector += latticelist[i]
        for j in range(len(newlist)):
            nextvector -= (np.dot(w,newlist[j])/(norm(newlist[j],1))) * newlist[j]
        newlist.append(nextvector)
    return newlist


def lovasz(ortholist):
    '''
    

    Parameters
    ----------
    ortholist : list 
        gram schmidts orthogonal list

    Returns
    -------
    ortholist : list
        return ortholist that satisfy lovasz condition
    boolean : bool
        if it is true, there is a swap
        otherwise, no swap

    '''
    boolean = False
    for k in range(len(ortholist)-1):
        if 0.75*norm(ortholist[k],1)> norm(ortholist[k+1]+(np.dot(ortholist[k+1],ortholist[k])/(norm(ortholist[k+1],1)))*ortholist[k],1):
            temp = ortholist[k]
            ortholist[k] = ortholist[k+1]
            ortholist[k+1] = temp 
            boolean = True
            break
    return ortholist, boolean


def lllalgo(l):
   if lovasz(gramschmidts(l))[1] == True:
       return lllalgo(lovasz(gramschmidts(l))[0])
   else:
       return lovasz(gramschmidts(l))[0]
        
print(lllalgo(latticelist))   


     

    