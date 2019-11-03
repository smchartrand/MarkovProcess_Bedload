# -*- coding: utf-8 -*-
"""
Created on Sun Jun 24 19:44:13 2018

@author: schartrand
"""
from numba import jit
import numpy as np


@jit(nopython=True)
def find_first(item, vec):
    """return the index of the first occurence of item in vec"""
    for i in xrange(len(vec)):
        if vec[i] > item:
            return i
    return -1

idx = 700  
a = np.arange(1, 1000, dtype=None)
c = find_first(idx,a)