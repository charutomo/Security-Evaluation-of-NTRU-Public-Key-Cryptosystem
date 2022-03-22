# -*- coding: utf-8 -*-
"""
Created on Mon Mar 21 13:59:02 2022

@author: Charissa
"""

import numpy as np
def convert(txtfile):
    '''
    

    Parameters
    ----------
    txtfile : str
        provide the name of text file with lattice in the form
        0,1,2      converts to list of basis
        3,4,5    -----------------------------> [np.array([0,1,2)]),np.array([3,4,5)]), np.array([6,7,8)])]
        6,7,8 

    Returns
    -------
    latticelist : list
        output list of arrays in its rows to perform gram schmidts orthogonalisation afterwards



    '''
    
    with open(str(txtfile), 'r') as f:
         lattice = np.array([[int(value) for value in line.split(',')] for line in f])
    
    latticelist = []
    for i in range(len(lattice)):
            latticelist.append(lattice[i])
    return latticelist
