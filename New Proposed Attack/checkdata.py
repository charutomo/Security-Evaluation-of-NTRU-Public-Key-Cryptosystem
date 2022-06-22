# -*- coding: utf-8 -*-
"""
Created on Fri Mar  4 19:02:13 2022

@author: Charissa
"""

import numpy as np
from numpy.linalg import norm 
from conversiontxtfile import converttolist

latlist = converttolist('answerdata.txt')
orig = converttolist("testdata.txt")


def SVFinder(latlist):
    val = 0
    lt =[]
    for i in range(len(latlist)):
        for j in range(len(latlist[i])):
            val += latlist[i][j]*orig[j]
        lt.append(val)
        val = 0
    normlt = sorted(lt, key = norm)
    return normlt[0], norm(normlt[0],2)
        

def SV(latlist):
    largestnorm =[]
    for i in range(len(latlist)):
        largestnorm.append(norm(latlist[i],2))
    return (min(largestnorm))

SV(latlist)
SVFinder(latlist)
