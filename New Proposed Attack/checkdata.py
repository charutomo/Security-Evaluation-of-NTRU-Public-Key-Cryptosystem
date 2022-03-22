# -*- coding: utf-8 -*-
"""
Created on Fri Mar  4 19:02:13 2022

@author: Charissa
"""

import numpy as np
from numpy.linalg import norm 
from conversiontxtfile import convert 

latlist = convert('answerdata.txt')
print(latlist)

def SV(latlist):
    largestnorm =[]
    for i in range(len(latlist)):
        largestnorm.append(norm(latlist[i],2))
    return (max(largestnorm))

print(SV(latlist))